#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
r'''! Post-solve coil limit checks for TokaMaker equilibria

Given a solved equilibrium and user-supplied engineering limits, report each
check with its value, its limit and the fraction of the limit used, and
identify which limit binds. No machine limits are included here: all limits
are supplied by the user.

Limits are given as either a scalar `L` (meaning `|value| <= L`) or a pair
`[min, max]` with `min <= 0 <= max` (either entry may be `None` for no limit).
A limit of `None` skips that check. The fraction of a limit used is
`value/max` for positive values and `value/min` for negative values, so a
value above 1 means the limit is exceeded.

Quantities and units:
 - `current`: current per turn [A] (the coil current used by TokaMaker,
   see @ref TokaMaker.TokaMaker.get_coil_currents "get_coil_currents")
 - `ampere_turns`: current times number of turns [A]
 - `current_density`: toroidal current density averaged over the coil
   cross-section [A/m^2]
 - `force_r`: net radial body force on the coil, integrated over the full torus [N]
   (\f$ 2 \pi \int J_\phi B_z R dA \f$; divide by \f$ 2 \pi \f$ for N/rad)
 - `force_z`: net vertical force on the coil [N] (\f$ -2 \pi \int J_\phi B_r R dA \f$)
 - `energy_coils`, `energy_total`: poloidal magnetic energy of the coils only,
   and of coils plus plasma, from the inductance matrix [J]

@authors Eleanor Adams
@date October 2026
@ingroup doxy_oft_python
'''
from warnings import warn
import numpy

## Quantities that can be limited for each coil or coil set
COIL_QUANTITIES = ('current', 'ampere_turns', 'current_density', 'force_r', 'force_z')
## Quantities that can be limited for a group of coils (weighted sums)
GROUP_QUANTITIES = ('ampere_turns', 'force_r', 'force_z')
## Quantities that can be limited for the whole coil system
GLOBAL_QUANTITIES = ('energy_coils', 'energy_total')
## Units of each quantity, for display
UNITS = {
    'current': 'A',
    'ampere_turns': 'A-turn',
    'current_density': 'A/m^2',
    'force_r': 'N',
    'force_z': 'N',
    'energy_coils': 'J',
    'energy_total': 'J',
}


def _get_bounds(limit):
    '''! Convert a limit specification to a (min, max) pair

    @param limit Scalar `L` (meaning `[-L, L]`) or pair `[min, max]` (`None` entries mean no limit)
    @result (min, max) with `min <= 0 <= max`
    '''
    if numpy.ndim(limit) == 0:
        try:
            value = float(limit)
        except (TypeError, ValueError):
            raise ValueError('Limits must be a number or a pair [min, max], got {0!r}'.format(limit))
        if not value > 0.0:
            raise ValueError('Scalar limits must be positive, got {0}'.format(limit))
        return (-value, value)
    if numpy.ndim(limit) != 1 or len(limit) != 2:
        raise ValueError('Limits must be a number or a pair [min, max], got {0!r}'.format(limit))
    try:
        lower = -numpy.inf if limit[0] is None else float(limit[0])
        upper = numpy.inf if limit[1] is None else float(limit[1])
    except (TypeError, ValueError):
        raise ValueError('Limit pair entries must be numbers or None, got {0!r}'.format(limit))
    if (not lower <= 0.0) or (not upper >= 0.0):
        raise ValueError('Limit pairs must satisfy min <= 0 <= max, got {0}'.format(limit))
    return (lower, upper)


def _fraction_used(value, bounds):
    '''! Compute the fraction of a limit used by a value

    @param value Value to check
    @param bounds (min, max) pair from @ref _get_bounds
    @result Fraction of the active limit used (> 1 means the limit is exceeded), active limit
    '''
    lower, upper = bounds
    if value > 0.0:
        if upper == 0.0:
            return numpy.inf, upper
        return value/upper, upper
    if value < 0.0:
        if lower == 0.0:
            return numpy.inf, lower
        return value/lower, lower
    return 0.0, upper


def _make_check(name, kind, quantity, value, limit):
    '''! Build a single check record

    @param name Name of coil, coil set, group, or "system"
    @param kind One of "coil", "group", or "global"
    @param quantity Quantity name (see @ref UNITS)
    @param value Computed value (must be finite)
    @param limit User limit (scalar or [min, max])
    @result Check record dictionary
    '''
    if not numpy.isfinite(value):
        raise ValueError('Value of {0} for "{1}" is not finite ({2}), cannot check against limit'.format(quantity, name, value))
    bounds = _get_bounds(limit)
    fraction, active_limit = _fraction_used(value, bounds)
    return {
        'name': name,
        'kind': kind,
        'quantity': quantity,
        'units': UNITS[quantity],
        'value': float(value),
        'bounds': bounds,
        'limit': active_limit,
        'fraction_used': fraction,
        'ok': bool(fraction <= 1.0),
    }


def _lookup_coil(quantities, name):
    '''! Get computed values for a coil set or a single coil region by name

    @param quantities Output of @ref compute_coil_quantities
    @param name Coil set name or coil (region) name
    @result Dictionary of values for that coil
    '''
    if name in quantities['coil_sets']:
        return quantities['coil_sets'][name]
    if name in quantities['coils']:
        return quantities['coils'][name]
    raise KeyError('Unknown coil "{0}"'.format(name))


def _warn_double_counting(quantities, group_name, names):
    '''! Warn if a group lists both a coil set and one of its regions

    @param quantities Output of @ref compute_coil_quantities
    @param group_name Name of group (for the warning message)
    @param names Coil set and coil region names in the group
    '''
    for name in names:
        if (name in quantities['coil_sets']) or (name not in quantities['coils']):
            continue
        set_name = quantities['coils'][name]['coil_set']
        if (set_name != name) and (set_name in names):
            warn('Group "{0}" lists coil set "{1}" and its region "{2}", which counts "{2}" twice'.format(group_name, set_name, name))


def evaluate_limits(quantities, coil_limits=None, group_limits=None, global_limits=None):
    '''! Compare computed coil quantities against user limits

    This is the numerical part of @ref check_coil_limits and does not need an
    equilibrium, so it can be reused with stored quantities. Limits of `None`
    are skipped.

    @param quantities Output of @ref compute_coil_quantities
    @param coil_limits Limits per coil set or coil region, `{name: {quantity: limit}}`
    with quantities from @ref COIL_QUANTITIES
    @param group_limits Limits on weighted sums over coils,
    `{group_name: {'coils': [names] or {name: weight}, quantity: limit}}` with
    quantities from @ref GROUP_QUANTITIES
    @param global_limits Limits on system quantities, `{quantity: limit}` with
    quantities from @ref GLOBAL_QUANTITIES
    @result List of check record dictionaries
    '''
    checks = []
    for name, limits in (coil_limits or {}).items():
        values = _lookup_coil(quantities, name)
        for quantity, limit in limits.items():
            if quantity not in COIL_QUANTITIES:
                raise ValueError('Unknown coil limit "{0}" for "{1}", must be one of {2}'.format(quantity, name, COIL_QUANTITIES))
            if limit is not None:
                checks.append(_make_check(name, 'coil', quantity, values[quantity], limit))
    for group_name, group in (group_limits or {}).items():
        if 'coils' not in group:
            raise ValueError('Group "{0}" must list its coils under the "coils" key'.format(group_name))
        weights = group['coils']
        if not isinstance(weights, dict):
            weights = {name: 1.0 for name in weights}
        _warn_double_counting(quantities, group_name, list(weights.keys()))
        for quantity, limit in group.items():
            if quantity == 'coils':
                continue
            if quantity not in GROUP_QUANTITIES:
                raise ValueError('Unknown group limit "{0}" for "{1}", must be one of {2}'.format(quantity, group_name, GROUP_QUANTITIES))
            if limit is None:
                continue
            value = 0.0
            for name, weight in weights.items():
                value += weight*_lookup_coil(quantities, name)[quantity]
            checks.append(_make_check(group_name, 'group', quantity, value, limit))
    for quantity, limit in (global_limits or {}).items():
        if quantity not in GLOBAL_QUANTITIES:
            raise ValueError('Unknown global limit "{0}", must be one of {1}'.format(quantity, GLOBAL_QUANTITIES))
        if limit is not None:
            checks.append(_make_check('system', 'global', quantity, quantities[quantity], limit))
    return checks


def _check_supported(tMaker_obj):
    '''! Raise an error for configurations the coil limit checks do not support

    @param tMaker_obj TokaMaker object
    '''
    if tMaker_obj.settings.mirror_mode:
        raise NotImplementedError('Coil limit checks are not supported in mirror mode')
    if len(tMaker_obj.dist_coils) > 0:
        names = [tMaker_obj.coil_set_names[coil_id] for coil_id in tMaker_obj.dist_coils]
        raise NotImplementedError('Coil limit checks assume uniform current density, but coil(s) {0} '
                                  'have a non-uniform current distribution'.format(names))


def _plasma_current(tMaker_obj):
    r'''! Get the plasma current carried by the present \f$ \psi \f$

    `get_globals` reports the plasma current from the profiles, which can be
    nonzero after a vacuum solve even though \f$ \psi \f$ carries no plasma current.
    The reported value is used only if the current from \f$ \Delta^* \psi \f$ in the
    plasma region agrees with it to within 50%.

    @param tMaker_obj TokaMaker object
    @result Plasma current [A] (0 if \f$ \psi \f$ has no plasma current)
    '''
    Ip = tMaker_obj.get_globals()[0]
    if Ip == 0.0:
        return 0.0
    psi = tMaker_obj.get_psi(normalized=False)
    Ip_field = tMaker_obj.compute_area_integral(tMaker_obj.calc_delstar_curr(psi), reg_mask=1)
    if abs(Ip_field - Ip) > 0.5*abs(Ip):
        return 0.0
    return Ip


def compute_coil_forces(tMaker_obj):
    r'''! Compute net radial and vertical force on each coil region

    Uses the same field evaluation and sign conventions as
    @ref TokaMaker.util.compute_forces_components "compute_forces_components"
    (force density \f$ f_r = J_\phi B_z \f$, \f$ f_z = -J_\phi B_r \f$), applied to
    coil regions with a uniform current density. Forces are integrated over the
    full torus, so `force_r` is the total outward radial body force [N]
    (divide by \f$ 2 \pi \f$ for N/rad) and `force_z` is the net vertical force [N].
    The field includes all sources (plasma, coils, and any conductor currents in \f$ \psi \f$).

    @param tMaker_obj TokaMaker object with a solved equilibrium
    @result `{coil_name: {'force_r': F_r, 'force_z': F_z}}` per coil region [N]
    '''
    _check_supported(tMaker_obj)
    _, region_currents = tMaker_obj.get_coil_currents()
    field_eval = tMaker_obj.get_field_eval('B')
    ones = numpy.ones((tMaker_obj.np,))
    forces = {}
    for name, coil in tMaker_obj._coil_dict.items():
        reg_id = coil['reg_id']
        # Evaluate B only at the vertices of this coil region
        pt_mask = numpy.zeros((tMaker_obj.np,), dtype=bool)
        pt_mask[tMaker_obj.lc[tMaker_obj.reg == reg_id, :]] = True
        B = numpy.zeros((tMaker_obj.np, 3))
        B[pt_mask, :] = field_eval.eval(tMaker_obj.r[pt_mask, :2])
        area = tMaker_obj.compute_area_integral(ones, reg_mask=reg_id)
        J_phi = region_currents[reg_id-1]/area
        R = tMaker_obj.r[:, 0]
        forces[name] = {
            'force_r': 2.0*numpy.pi*J_phi*tMaker_obj.compute_area_integral(R*B[:, 2], reg_mask=reg_id),
            'force_z': -2.0*numpy.pi*J_phi*tMaker_obj.compute_area_integral(R*B[:, 0], reg_mask=reg_id),
        }
    return forces


def compute_coil_quantities(tMaker_obj):
    r'''! Compute all coil quantities used by the limit checks

    Values are computed for each coil region (keyed by coil name, as in the mesh
    coil dictionary) and for each coil set (keyed by coil set name, as used by
    @ref TokaMaker.TokaMaker.set_coil_bounds "set_coil_bounds"). For a coil set,
    `ampere_turns` and forces are net sums over its regions and `current_density`
    is the region value with the largest magnitude.

    With no plasma current in \f$ \psi \f$ (eg. after a vacuum solve), `energy_total`
    equals `energy_coils`.

    @note Coils with non-uniform current distributions (see
    @ref TokaMaker.TokaMaker.set_coil_current_dist "set_coil_current_dist") and
    mirror mode are not supported and raise `NotImplementedError`.

    @param tMaker_obj TokaMaker object with a solved equilibrium
    @result Dictionary with keys `coils`, `coil_sets` (per-coil values, see @ref COIL_QUANTITIES),
    `energy_coils` and `energy_total` [J]
    '''
    _check_supported(tMaker_obj)
    set_currents, region_currents = tMaker_obj.get_coil_currents()
    forces = compute_coil_forces(tMaker_obj)
    ones = numpy.ones((tMaker_obj.np,))
    coils = {}
    coil_sets = {}
    for set_name, coil_set in tMaker_obj.coil_sets.items():
        current = set_currents[set_name]
        set_values = {
            'current': current,
            'ampere_turns': 0.0,
            'current_density': 0.0,
            'force_r': 0.0,
            'force_z': 0.0,
            'area': 0.0,
        }
        for sub_coil in coil_set['sub_coils']:
            reg_id = sub_coil['reg_id']
            area = tMaker_obj.compute_area_integral(ones, reg_mask=reg_id)
            ampere_turns = region_currents[reg_id-1]
            values = {
                'coil_set': set_name,
                'current': current,
                'nturns': sub_coil.get('nturns', 1.0),
                'ampere_turns': ampere_turns,
                'current_density': ampere_turns/area,
                'force_r': forces[sub_coil['name']]['force_r'],
                'force_z': forces[sub_coil['name']]['force_z'],
                'area': area,
            }
            coils[sub_coil['name']] = values
            set_values['ampere_turns'] += ampere_turns
            set_values['force_r'] += values['force_r']
            set_values['force_z'] += values['force_z']
            set_values['area'] += area
            if abs(values['current_density']) > abs(set_values['current_density']):
                set_values['current_density'] = values['current_density']
        coil_sets[set_name] = set_values
    # Stored energy from the inductance matrix (coil currents per turn, plus plasma current)
    current_vec = numpy.zeros((tMaker_obj.ncoils,))
    for set_name, coil_set in tMaker_obj.coil_sets.items():
        current_vec[coil_set['id']] = set_currents[set_name]
    energy_coils = 0.5*numpy.dot(current_vec, numpy.dot(tMaker_obj.Lcoils, current_vec))
    energy_total = energy_coils
    Ip = _plasma_current(tMaker_obj)
    if Ip != 0.0:
        # Plasma row of the inductance matrix is normalized by Ip, so only use it with plasma current
        Lmat = tMaker_obj.get_coil_Lmat()
        full_vec = numpy.append(current_vec, Ip)
        energy_total = 0.5*numpy.dot(full_vec, numpy.dot(Lmat, full_vec))
    return {
        'coils': coils,
        'coil_sets': coil_sets,
        'energy_coils': energy_coils,
        'energy_total': energy_total,
    }


class CoilLimitReport():
    '''! Results of a coil limit check (see @ref check_coil_limits)'''
    def __init__(self, checks, quantities=None):
        '''! Create report from check records

        @param checks List of check records from @ref evaluate_limits
        @param quantities Computed quantities from @ref compute_coil_quantities
        '''
        ## All limit checks, each a dictionary with `name`, `kind`, `quantity`, `units`,
        ## `value`, `bounds`, `limit` (active bound), `fraction_used`, and `ok`
        self.checks = checks
        ## Computed quantities used for the checks
        self.quantities = quantities
        ## Check with the largest fraction of its limit used (`None` if no checks)
        self.binding = None
        if len(checks) > 0:
            self.binding = max(checks, key=lambda check: check['fraction_used'])

    @property
    def ok(self):
        '''! True if all limit checks pass'''
        return all(check['ok'] for check in self.checks)

    @property
    def violations(self):
        '''! List of limit checks that exceed their limit'''
        return [check for check in self.checks if not check['ok']]

    def max_fraction_by_name(self):
        '''! Largest fraction of limit used for each coil, group, or "system"

        @result Dictionary `{name: (fraction_used, quantity)}`, in order of first appearance
        '''
        result = {}
        for check in self.checks:
            name = check['name']
            if (name not in result) or (check['fraction_used'] > result[name][0]):
                result[name] = (check['fraction_used'], check['quantity'])
        return result

    def summary(self):
        '''! Readable text table of all checks

        @result Summary text
        '''
        lines = []
        if self.binding is not None:
            lines.append('Binding limit: {0} {1} ({2:.1f}% of limit)'.format(
                self.binding['name'], self.binding['quantity'], 100.0*self.binding['fraction_used']))
        lines.append('All limits satisfied' if self.ok else 'LIMITS EXCEEDED')
        header = '{0:<16} {1:<16} {2:>12} {3:>12} {4:<8} {5:>8}'.format('Name', 'Quantity', 'Value', 'Limit', 'Units', 'Used')
        lines.append(header)
        lines.append('-'*len(header))
        for check in self.checks:
            flag = ''
            if not check['ok']:
                flag = ' !'
            elif check is self.binding:
                flag = ' *'
            lines.append('{0:<16} {1:<16} {2:>12.4e} {3:>12.4e} {4:<8} {5:>7.1f}%{6}'.format(
                check['name'], check['quantity'], check['value'], check['limit'],
                check['units'], 100.0*check['fraction_used'], flag))
        return '\n'.join(lines)

    def plot_margins(self, fig, ax):
        '''! Bar plot of the largest fraction of limit used for each coil, group, or "system"

        @param fig Figure to add to (unused)
        @param ax Axis to add to
        '''
        max_fractions = self.max_fraction_by_name()
        names = list(max_fractions.keys())
        fractions = numpy.array([max_fractions[name][0] for name in names])
        colors = ['tab:red' if fraction > 1.0 else 'tab:blue' for fraction in fractions]
        ax.bar(names, 100.0*fractions, color=colors)
        ax.axhline(100.0, color='k', linestyle='--')
        ax.set_ylabel('Largest fraction of limit used [%]')
        ax.tick_params(axis='x', labelrotation=90)


def check_coil_limits(tMaker_obj, coil_limits=None, group_limits=None, global_limits=None):
    '''! Check a solved equilibrium against user-supplied coil limits

    This is a post-solve check and is separate from the solver constraint set by
    @ref TokaMaker.TokaMaker.set_coil_bounds "set_coil_bounds", but uses the same
    coil (set) names. A coil set name refers to all regions in that set; a coil
    region name refers to that region only.

    Example:
    ~~~~~~~~~~{.py}
    report = check_coil_limits(mygs,
        coil_limits={'PF1': {'current': 45.E3, 'force_z': [-50.E6, 50.E6]}},
        group_limits={'CS_upper': {'coils': ['CS1U', 'CS2U'], 'force_z': 100.E6}},
        global_limits={'energy_coils': 5.E9})
    print(report.summary())
    ~~~~~~~~~~

    @param tMaker_obj TokaMaker object with a solved equilibrium
    @param coil_limits Limits per coil set or coil region, `{name: {quantity: limit}}`
    with quantities from @ref COIL_QUANTITIES
    @param group_limits Limits on weighted sums over coils,
    `{group_name: {'coils': [names] or {name: weight}, quantity: limit}}` with
    quantities from @ref GROUP_QUANTITIES
    @param global_limits Limits on system quantities, `{quantity: limit}` with
    quantities from @ref GLOBAL_QUANTITIES
    @result @ref CoilLimitReport
    '''
    quantities = compute_coil_quantities(tMaker_obj)
    checks = evaluate_limits(quantities, coil_limits, group_limits, global_limits)
    return CoilLimitReport(checks, quantities)
