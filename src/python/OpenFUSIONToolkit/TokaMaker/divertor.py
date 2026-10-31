#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
r'''! Divertor geometry metrics for TokaMaker equilibria

Computes strike points, flux expansion, connection lengths, target incidence angles,
\f$ \Delta R_{sep} \f$ and secondary X-point distances for a solved equilibrium and its
limiter (wall) contour. See @ref OpenFUSIONToolkit.TokaMaker.divertor.compute_divertor_metrics "compute_divertor_metrics"
for definitions and units.

@authors GitLostintheSauce
@date October 2026
@ingroup doxy_oft_python
'''
from warnings import warn
import numpy


class _TokaMakerField():
    r'''! Thin adapter exposing \f$ \psi \f$ and \f$ B \f$ of a TokaMaker equilibrium as vectorized functions'''
    def __init__(self, equil):
        '''! Create field adapter

        @param equil TokaMaker equilibrium object (@ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker_equilibrium "TokaMaker_equilibrium")
        '''
        self._equil = equil
        self._psi_eval = equil.get_field_eval('PSI')
        self._B_eval = equil.get_field_eval('B')

    def psi(self, pts):
        r'''! Evaluate poloidal flux

        @param pts Locations [n,2] (R, Z) [m]
        @result \f$ \psi \f$ [n] [Wb/rad] (NaN outside the mesh)
        '''
        return self._psi_eval.eval(numpy.atleast_2d(pts))[:, 0]

    def B(self, pts):
        '''! Evaluate magnetic field

        @param pts Locations [n,2] (R, Z) [m]
        @result (B_R, B_phi, B_Z) [n,3] [T] (NaN outside the mesh)
        '''
        return self._B_eval.eval(numpy.atleast_2d(pts))

    def psi_norm(self, psi):
        r'''! Convert \f$ \psi \f$ to normalized flux \f$ \hat{\psi} \f$ (0 on axis, 1 on the LCFS)

        @param psi Poloidal flux [Wb/rad]
        @result \f$ \hat{\psi} \f$
        '''
        return self._equil.abspsi_to_normalized(psi)


def _cross2(a, b):
    '''! z-component of the cross product of 2D vectors (broadcasts)'''
    return a[..., 0]*b[..., 1] - a[..., 1]*b[..., 0]


def _wall_segments(walls):
    '''! Stack all wall contours into segment arrays

    @param walls List of closed contours [n,2]
    @result Segment start points [m,2], segment vectors [m,2]
    '''
    starts = []
    vecs = []
    for contour in walls:
        starts.append(contour[:-1, :])
        vecs.append(numpy.diff(contour, axis=0))
    return numpy.vstack(starts), numpy.vstack(vecs)


def _first_wall_hit(p, q, segs):
    '''! Find the first intersection of segment p->q with the wall

    @param p Segment start [2]
    @param q Segment end [2]
    @param segs Wall segments from `_wall_segments`
    @result Fraction along p->q of the first hit and the wall segment index, or (None, -1)
    '''
    starts, vecs = segs
    r = q - p
    denom = _cross2(r, vecs)
    w = starts - p
    with numpy.errstate(divide='ignore', invalid='ignore'):
        s = _cross2(w, vecs)/denom
        u = _cross2(w, r)/denom
    valid = (numpy.abs(denom) > 1.E-14) & (s > 1.E-12) & (s <= 1.0) & (u >= 0.0) & (u <= 1.0)
    if not valid.any():
        return None, -1
    idx = numpy.flatnonzero(valid)
    i = idx[s[idx].argmin()]
    return s[i], i


def _point_in_contour(pt, contour):
    '''! Ray-casting point-in-polygon test

    @param pt Point [2]
    @param contour Closed contour [n,2]
    @result True if `pt` is inside `contour`
    '''
    x0, y0 = contour[:-1, 0], contour[:-1, 1]
    x1, y1 = contour[1:, 0], contour[1:, 1]
    crosses = (y0 > pt[1]) != (y1 > pt[1])
    with numpy.errstate(divide='ignore', invalid='ignore'):
        x_int = x0 + (pt[1] - y0)*(x1 - x0)/(y1 - y0)
    return bool(numpy.count_nonzero(crosses & (pt[0] < x_int)) % 2)


def _distance_to_segments(pt, segs):
    '''! Minimum distance from a point to a set of segments

    @param pt Point [2]
    @param segs Segments from `_wall_segments`
    @result Distance [m]
    '''
    starts, vecs = segs
    len2 = numpy.maximum((vecs**2).sum(axis=1), 1.E-30)
    t = numpy.clip(((pt - starts)*vecs).sum(axis=1)/len2, 0.0, 1.0)
    closest = starts + t[:, None]*vecs
    return numpy.sqrt(((closest - pt)**2).sum(axis=1)).min()


def _bisect(fun, a, b, fa, niter=60):
    '''! Bisection for a sign change of `fun` on [a, b]

    @param fun Scalar function
    @param a Left end
    @param b Right end
    @param fa fun(a)
    @param niter Number of iterations
    @result Root location
    '''
    for _ in range(niter):
        c = 0.5*(a + b)
        fc = fun(c)
        if fc == 0.0:
            return c
        if (fc > 0.0) == (fa > 0.0):
            a, fa = c, fc
        else:
            b = c
    return 0.5*(a + b)


def _omp_crossing(field, psi_target, o_point, R_wall, nsample=400):
    r'''! Find where a flux surface crosses the outboard midplane

    @param field Field adapter
    @param psi_target \f$ \psi \f$ of the surface [Wb/rad]
    @param o_point Magnetic axis (R, Z) [m]
    @param R_wall Radius where the midplane ray hits the wall [m]
    @param nsample Number of samples used to bracket the crossing
    @result Innermost crossing radius [m] (NaN if the surface does not reach the midplane)
    '''
    R = numpy.linspace(o_point[0] + 0.01*(R_wall - o_point[0]), R_wall, nsample)
    pts = numpy.column_stack((R, numpy.full(nsample, o_point[1])))
    f = field.psi(pts) - psi_target
    for i in range(nsample - 1):
        if not (numpy.isfinite(f[i]) and numpy.isfinite(f[i+1])):
            continue
        if f[i] == 0.0:
            return R[i]
        if f[i]*f[i+1] < 0.0:
            def fun(r):
                return field.psi(numpy.r_[r, o_point[1]])[0] - psi_target
            return _bisect(fun, R[i], R[i+1], f[i])
    return numpy.nan


def _separatrix_branches(field, x_point, psi_x, radius, nsample=360):
    r'''! Find separatrix branch directions around an X-point

    Samples \f$ \psi - \psi_X \f$ on a circle around the X-point and locates its zeros.

    @param field Field adapter
    @param x_point X-point (R, Z) [m]
    @param psi_x \f$ \psi \f$ at the X-point [Wb/rad]
    @param radius Circle radius [m]
    @param nsample Number of samples on the circle
    @result Unit direction vectors of each branch [nbranch,2]
    '''
    theta = numpy.linspace(0.0, 2.0*numpy.pi, nsample, endpoint=False)
    dtheta = theta[1]
    pts = x_point + radius*numpy.column_stack((numpy.cos(theta), numpy.sin(theta)))
    f = field.psi(pts) - psi_x

    def fun(th):
        return field.psi(x_point + radius*numpy.r_[numpy.cos(th), numpy.sin(th)])[0] - psi_x
    dirs = []
    for i in range(nsample):
        j = (i + 1) % nsample
        if not (numpy.isfinite(f[i]) and numpy.isfinite(f[j])):
            continue
        if f[i]*f[j] < 0.0:
            th = _bisect(fun, theta[i], theta[i] + dtheta, f[i], niter=40)
            dirs.append([numpy.cos(th), numpy.sin(th)])
    return numpy.array(dirs).reshape((-1, 2))


def _trace_field_line(field, start, psi_target, direction, segs, x_points, step, max_length, stop_radius=0.0):
    r'''! Follow a field line in the poloidal plane until it reaches the wall

    Integrates \f$ d\mathbf{x}/dl_{pol} = \pm \mathbf{B}_p/B_p \f$ and \f$ dL_\parallel/dl_{pol} = |B|/B_p \f$
    with RK4, projecting back onto the target flux surface after each step. The step is reduced
    near X-points, where the field direction changes rapidly.

    @param field Field adapter
    @param start Starting point (R, Z) [m]
    @param psi_target \f$ \psi \f$ of the surface being followed [Wb/rad]
    @param direction +1 to follow \f$ +\mathbf{B}_p \f$, -1 to follow \f$ -\mathbf{B}_p \f$
    @param segs Wall segments from `_wall_segments`
    @param x_points X-points used for step control and stopping [n,2]
    @param step Maximum poloidal step [m]
    @param max_length Maximum poloidal length to trace [m]
    @param stop_radius Stop if the line comes within this distance of an X-point (0 to disable) [m]
    @result Dictionary with `pts` [n,2], cumulative `L_pol` [n] and `L_par` [n] [m], `status`
    ('wall', 'xpoint', 'mesh' or 'max_length') and `wall_segment` index of the hit wall segment (-1 if none)
    '''
    def rhs(x):
        B = field.B(x)[0]
        Bp = numpy.hypot(B[0], B[2])
        if (not numpy.all(numpy.isfinite(B))) or (Bp <= 0.0):
            return None
        return numpy.r_[direction*B[0]/Bp, direction*B[2]/Bp, numpy.linalg.norm(B)/Bp]

    def xpt_dist(x):
        if x_points.shape[0] == 0:
            return numpy.inf
        return numpy.sqrt(((x_points - x)**2).sum(axis=1)).min()

    x = numpy.array(start, dtype=numpy.float64)
    pts = [x.copy()]
    L_pol = [0.0]
    L_par = [0.0]
    status = 'max_length'
    wall_seg = -1
    k1 = rhs(x)
    if k1 is None:
        status = 'mesh'
    while (k1 is not None) and (L_pol[-1] < max_length):
        h = min(step, max(0.02*step, 0.25*xpt_dist(x)))
        k2 = rhs(x + 0.5*h*k1[:2])
        k3 = rhs(x + 0.5*h*k2[:2]) if k2 is not None else None
        k4 = rhs(x + h*k3[:2]) if k3 is not None else None
        if k4 is None:
            status = 'mesh'
            break
        incr = h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        x_new = x + incr[:2]
        # Project back onto the target surface to remove integration drift
        psi_new = field.psi(x_new)[0]
        B_new = field.B(x_new)[0]
        grad = x_new[0]*numpy.r_[B_new[2], -B_new[0]]
        grad2 = numpy.dot(grad, grad)
        if numpy.isfinite(psi_new) and (grad2 > 0.0):
            corr = (psi_new - psi_target)*grad/grad2
            if numpy.linalg.norm(corr) < 0.5*h:
                x_new = x_new - corr
        # Check for wall intersection
        s, iseg = _first_wall_hit(x, x_new, segs)
        if s is not None:
            pts.append(x + s*(x_new - x))
            L_pol.append(L_pol[-1] + s*h)
            L_par.append(L_par[-1] + s*incr[2])
            status = 'wall'
            wall_seg = iseg
            break
        pts.append(x_new)
        L_pol.append(L_pol[-1] + h)
        L_par.append(L_par[-1] + incr[2])
        x = x_new
        if (stop_radius > 0.0) and (L_pol[-1] > 4.0*stop_radius) and (xpt_dist(x) < stop_radius):
            status = 'xpoint'
            break
        k1 = rhs(x)
        if k1 is None:
            status = 'mesh'
    return {
        'pts': numpy.array(pts),
        'L_pol': numpy.array(L_pol),
        'L_par': numpy.array(L_par),
        'status': status,
        'wall_segment': wall_seg
    }


def _polyline_distance(pt, pts):
    '''! Minimum distance from a point to the vertices of a densely sampled polyline'''
    return numpy.sqrt(((pts - pt)**2).sum(axis=1)).min()


def _divertor_metrics(field, o_point, x_points, walls, sol_width=3.E-3, step=5.E-3, max_length=None, diverted=None):
    r'''! Compute divertor geometry metrics from a field adapter

    Core of @ref OpenFUSIONToolkit.TokaMaker.divertor.compute_divertor_metrics "compute_divertor_metrics", independent of the TokaMaker object so it can be used with any
    field providing `psi(pts)`, `B(pts)` and optionally `psi_norm(psi)`.

    @param field Field adapter
    @param o_point Magnetic axis (R, Z) [m]
    @param x_points Candidate X-points [n,2] [m]
    @param walls List of wall contours [n,2] [m]
    @param sol_width Outboard-midplane distance outside the separatrix of the SOL surface used for connection lengths [m]
    @param step Maximum poloidal step for field-line tracing [m]
    @param max_length Maximum poloidal length of a single trace [m] (default: 10x total wall length)
    @param diverted Diverted flag reported in the output (informational)
    @result Metrics dictionary (see @ref OpenFUSIONToolkit.TokaMaker.divertor.compute_divertor_metrics "compute_divertor_metrics")
    '''
    o_point = numpy.asarray(o_point, dtype=numpy.float64)
    closed_walls = []
    for contour in walls:
        contour = numpy.asarray(contour, dtype=numpy.float64)
        if numpy.linalg.norm(contour[0, :] - contour[-1, :]) > 0.0:
            contour = numpy.vstack((contour, contour[:1, :]))
        closed_walls.append(contour)
    segs = _wall_segments(closed_walls)
    if max_length is None:
        max_length = 10.0*numpy.linalg.norm(segs[1], axis=1).sum()
    if x_points is None:
        x_points = numpy.zeros((0, 2))
    x_points = numpy.asarray(x_points, dtype=numpy.float64).reshape((-1, 2))
    has_norm = hasattr(field, 'psi_norm')

    # Outboard midplane ray
    s, _ = _first_wall_hit(o_point, o_point + numpy.r_[1.E3, 0.0], segs)
    if s is None:
        raise ValueError('Outboard midplane ray from the magnetic axis does not intersect the wall')
    R_wall_omp = o_point[0] + s*1.E3

    # Characterize X-points
    psi_x = field.psi(x_points) if x_points.shape[0] > 0 else numpy.zeros((0,))
    xpt_info = []
    for i in range(x_points.shape[0]):
        inside = any(_point_in_contour(x_points[i, :], contour) for contour in closed_walls)
        R_omp = _omp_crossing(field, psi_x[i], o_point, R_wall_omp) if numpy.isfinite(psi_x[i]) else numpy.nan
        xpt_info.append({
            'R': x_points[i, 0],
            'Z': x_points[i, 1],
            'psi': psi_x[i],
            'psi_n': field.psi_norm(psi_x[i]) if has_norm else numpy.nan,
            'inside_wall': inside,
            'R_OMP': R_omp,
            'distance_to_wall': _distance_to_segments(x_points[i, :], segs),
            'primary': False
        })
    inside_idx = [i for i, info in enumerate(xpt_info) if info['inside_wall']]
    candidates = [i for i in inside_idx if numpy.isfinite(xpt_info[i]['R_OMP'])]
    if len(candidates) == 0:
        raise ValueError('No X-point inside the wall whose separatrix reaches the outboard midplane')
    ipri = min(candidates, key=lambda i: xpt_info[i]['R_OMP'])
    xpt_info[ipri]['primary'] = True
    x_inside = x_points[inside_idx, :]

    # Outboard midplane quantities on the primary separatrix
    R_omp = xpt_info[ipri]['R_OMP']
    B_omp = field.B(numpy.r_[R_omp, o_point[1]])[0]

    # Double-null balance: lower minus upper separatrix radius at the OMP
    upper = xpt_info[ipri]['Z'] > o_point[1]
    opposite = [i for i, info in enumerate(xpt_info)
                if ((info['Z'] > o_point[1]) != upper) and numpy.isfinite(info['R_OMP'])]
    dR_sep = numpy.nan
    if len(opposite) > 0:
        iopp = min(opposite, key=lambda i: xpt_info[i]['R_OMP'])
        if upper:
            dR_sep = xpt_info[iopp]['R_OMP'] - R_omp
        else:
            dR_sep = R_omp - xpt_info[iopp]['R_OMP']

    # Trace separatrix legs from each X-point inside the wall (primary first)
    targets = {}
    for ix in [ipri] + [i for i in inside_idx if i != ipri]:
        x_pt = x_points[ix, :]
        others = numpy.array([x_points[j, :] for j in inside_idx if j != ix]).reshape((-1, 2))
        radius = 0.02
        if others.shape[0] > 0:
            radius = min(radius, 0.25*numpy.sqrt(((others - x_pt)**2).sum(axis=1)).min())
        radius = min(radius, 0.5*xpt_info[ix]['distance_to_wall'])
        branches = _separatrix_branches(field, x_pt, xpt_info[ix]['psi'], radius)
        if branches.shape[0] != 4:
            warn('Found {0} separatrix branches around X-point {1} (expected 4), skipping its legs'.format(branches.shape[0], ix))
            continue
        # Divertor legs point away from the magnetic axis; the other two branches wrap the core
        away = (x_pt - o_point)/numpy.linalg.norm(x_pt - o_point)
        legs = branches[numpy.argsort(branches.dot(away))[-2:], :]
        side = 'upper' if x_pt[1] > o_point[1] else 'lower'
        R_u = xpt_info[ix]['R_OMP']
        B_u = field.B(numpy.r_[R_u, o_point[1]])[0] if numpy.isfinite(R_u) else numpy.full(3, numpy.nan)
        Bp_u = numpy.hypot(B_u[0], B_u[2])
        for k in range(2):
            leg_dir = legs[k, :]
            name = '{0}_{1}'.format(side, 'outer' if leg_dir[0] >= legs[1-k, 0] else 'inner')
            if name in targets:
                name = '{0}_x{1}'.format(name, ix)
            start = x_pt + radius*leg_dir
            B_start = field.B(start)[0]
            direction = 1.0 if (B_start[0]*leg_dir[0] + B_start[2]*leg_dir[1]) >= 0.0 else -1.0
            trace = _trace_field_line(field, start, xpt_info[ix]['psi'], direction, segs,
                                      x_inside, step, max_length, stop_radius=0.5*radius)
            target = {
                'xpoint': ix,
                'primary': ix == ipri,
                'status': trace['status'],
                'trace': trace['pts'],
                'R_target': numpy.nan,
                'Z_target': numpy.nan,
                'leg_length_pol': numpy.nan,
                'B_pol_target': numpy.nan,
                'B_target': numpy.nan,
                'R_OMP': R_u,
                'poloidal_flux_expansion': numpy.nan,
                'toroidal_flux_expansion': numpy.nan,
                'total_flux_expansion': numpy.nan,
                'angle_poloidal_target': numpy.nan,
                'angle_of_incidence_target': numpy.nan,
                'xpoint_distance': numpy.nan,
                'secondary_xpoint': -1,
                'secondary_xpoint_distance': numpy.nan,
                'connection_length_target': numpy.nan,
                'connection_length_divertor': numpy.nan,
                'poloidal_connection_length_target': numpy.nan
            }
            targets[name] = target
            if trace['status'] != 'wall':
                continue
            pt = trace['pts'][-1, :]
            # Evaluate the field just inside the wall to avoid boundary round-off
            B_t = field.B(pt + 1.E-6*(trace['pts'][-2, :] - pt))[0]
            Bp_t = numpy.hypot(B_t[0], B_t[2])
            wall_vec = segs[1][trace['wall_segment'], :]
            normal = numpy.r_[wall_vec[1], -wall_vec[0]]/numpy.linalg.norm(wall_vec)
            Bn = abs(B_t[0]*normal[0] + B_t[2]*normal[1])
            target['R_target'] = pt[0]
            target['Z_target'] = pt[1]
            target['leg_length_pol'] = radius + trace['L_pol'][-1]
            target['B_pol_target'] = Bp_t
            target['B_target'] = numpy.linalg.norm(B_t)
            target['poloidal_flux_expansion'] = (R_u*Bp_u)/(pt[0]*Bp_t)
            target['toroidal_flux_expansion'] = pt[0]/R_u
            target['total_flux_expansion'] = Bp_u/Bp_t
            target['angle_poloidal_target'] = numpy.degrees(numpy.arcsin(min(1.0, Bn/Bp_t)))
            target['angle_of_incidence_target'] = numpy.degrees(numpy.arcsin(min(1.0, Bn/numpy.linalg.norm(B_t))))
            target['xpoint_distance'] = numpy.linalg.norm(pt - x_pt)
            for j in inside_idx:
                if j == ix:
                    continue
                dist = numpy.linalg.norm(pt - x_points[j, :])
                if not (dist >= target['secondary_xpoint_distance']):
                    target['secondary_xpoint'] = j
                    target['secondary_xpoint_distance'] = dist

    # Trace the SOL surface from the outboard midplane in both directions
    R_sol = R_omp + sol_width
    start = numpy.r_[R_sol, o_point[1]]
    psi_sol = field.psi(start)[0]
    sol = {
        'width': sol_width,
        'R_OMP': R_sol,
        'psi': psi_sol,
        'psi_n': field.psi_norm(psi_sol) if has_norm else numpy.nan,
        'ends': []
    }
    assign_tol = max(0.1, 30.0*sol_width)
    for direction in (1.0, -1.0):
        trace = _trace_field_line(field, start, psi_sol, direction, segs, x_inside, step, max_length)
        end = {
            'direction': direction,
            'status': trace['status'],
            'trace': trace['pts'],
            'L_par': trace['L_par'],
            'L_pol': trace['L_pol'],
            'R_target': numpy.nan,
            'Z_target': numpy.nan,
            'target': None,
            'connection_length': numpy.nan,
            'poloidal_length': numpy.nan
        }
        sol['ends'].append(end)
        if trace['status'] != 'wall':
            continue
        pt = trace['pts'][-1, :]
        end['R_target'] = pt[0]
        end['Z_target'] = pt[1]
        end['connection_length'] = trace['L_par'][-1]
        end['poloidal_length'] = trace['L_pol'][-1]
        # Assign to the nearest separatrix leg reaching the wall
        best_name = None
        best_dist = assign_tol
        for name, target in targets.items():
            if target['status'] != 'wall':
                continue
            dist = _polyline_distance(pt, target['trace'])
            if dist < best_dist:
                best_name, best_dist = name, dist
        if best_name is None:
            continue
        end['target'] = best_name
        target = targets[best_name]
        x_pt = x_points[target['xpoint'], :]
        ixpt = numpy.sqrt(((trace['pts'] - x_pt)**2).sum(axis=1)).argmin()
        target['connection_length_target'] = trace['L_par'][-1]
        target['poloidal_connection_length_target'] = trace['L_pol'][-1]
        target['connection_length_divertor'] = trace['L_par'][-1] - trace['L_par'][ixpt]

    return {
        'diverted': diverted,
        'R_axis': o_point[0],
        'Z_axis': o_point[1],
        'R_OMP': R_omp,
        'Z_OMP': o_point[1],
        'B_pol_OMP': numpy.hypot(B_omp[0], B_omp[2]),
        'B_tor_OMP': B_omp[1],
        'B_OMP': numpy.linalg.norm(B_omp),
        'psi_sep': xpt_info[ipri]['psi'],
        'dR_sep': dR_sep,
        'primary_xpoint': ipri,
        'xpoints': xpt_info,
        'targets': targets,
        'sol': sol,
        'wall': closed_walls
    }


def compute_divertor_metrics(tokamaker, sol_width=3.E-3, wall=None, step=5.E-3, max_length=None):
    r'''! Compute divertor geometry metrics for a solved TokaMaker equilibrium

    All lengths are in [m], fields in [T], \f$ \psi \f$ in [Wb/rad] and angles in [degrees].
    "OMP" denotes the outboard midplane, the horizontal ray \f$ Z = Z_{axis} \f$, \f$ R > R_{axis} \f$.

    **Algorithm**
    - The primary X-point is the X-point inside the wall whose separatrix crosses the OMP at the smallest
      radius (the active X-point when diverted). Secondary X-points inside the wall are also analyzed.
    - Separatrix branches are located on a small circle around each X-point. The two branches pointing away
      from the magnetic axis are the divertor legs; they are followed to the wall to give the strike points.
      Legs are named `{upper,lower}_{inner,outer}` (outer = the leg leaving the X-point at larger R).
      Names for secondary X-points that collide with an existing name get an `_x{index}` suffix.
    - Connection lengths use the SOL flux surface `sol_width` outside the primary separatrix at the OMP,
      traced in both directions from the OMP to the wall. Each end is attributed to the separatrix leg
      whose trace passes nearest to it (within max(0.1 m, 30*sol_width)).

    **Definitions** (u = OMP point of the same flux surface, t = target/strike point)
    - `poloidal_flux_expansion` \f$ f_{x} = (R_u B_{p,u})/(R_t B_{p,t}) \f$: ratio of flux-surface spacing at the
      target to that at the OMP. With \f$ F = R B_\phi \f$ constant in the SOL this equals the pitch-angle
      ratio \f$ (B_p/B_\phi)_u / (B_p/B_\phi)_t \f$.
    - `toroidal_flux_expansion` \f$ R_t/R_u \f$ (TORAX extended Lengyel `toroidal_flux_expansion`).
    - `total_flux_expansion` \f$ f_x R_t/R_u = B_{p,u}/B_{p,t} \f$.
    - `angle_poloidal_target`: angle between the separatrix leg and the wall in the poloidal plane.
    - `angle_of_incidence_target`: angle between the total field and the wall surface,
      \f$ \sin^{-1}(|\mathbf{B} \cdot \hat{n}|/|B|) \f$ (TORAX `angle_of_incidence_target`).
    - `connection_length_target`: parallel length \f$ L_\parallel \f$ from the OMP to the target along the SOL surface
      (TORAX `connection_length_target`).
    - `connection_length_divertor`: \f$ L_\parallel \f$ from the point of the SOL surface closest to the leg's
      X-point to the target (TORAX `connection_length_divertor`).
    - `dR_sep` \f$ \Delta R_{sep} = R_{OMP}(\psi_{X,lower}) - R_{OMP}(\psi_{X,upper}) \f$ using the primary X-point and
      the innermost X-point on the opposite side of the midplane; positive when the upper X-point is primary
      (DIII-D convention), NaN if there is no opposite X-point.
    - Flux expansions and angles are evaluated on the separatrix (\f$ \psi = \psi_X \f$) of each leg.

    **Returned dictionary**
    - `diverted`, `R_axis`, `Z_axis`, `Z_OMP`
    - `R_OMP`, `B_pol_OMP`, `B_tor_OMP`, `B_OMP`: primary separatrix at the OMP
    - `psi_sep`, `dR_sep`, `primary_xpoint` (index into `xpoints`)
    - `xpoints`: list of dicts with `R`, `Z`, `psi`, `psi_n`, `inside_wall`, `primary`, `R_OMP` (OMP radius of its
      separatrix, NaN if it does not reach the OMP), `distance_to_wall`
    - `targets`: dict keyed by leg name, each with `xpoint`, `primary`, `status` ('wall' if a strike point was
      found), `R_target`, `Z_target`, `leg_length_pol` (poloidal X-point to target length), `B_pol_target`,
      `B_target`, `R_OMP`, the flux expansions and angles above, `xpoint_distance` (strike point to own X-point),
      `secondary_xpoint`/`secondary_xpoint_distance` (nearest other X-point inside the wall to the strike point;
      -1/NaN if none), `connection_length_target`, `connection_length_divertor`,
      `poloidal_connection_length_target` (NaN if no SOL end reaches the target) and `trace` [n,2]
    - `sol`: dict with `width`, `R_OMP`, `psi`, `psi_n` and `ends` (one per tracing direction, with `target`,
      `R_target`, `Z_target`, `connection_length`, `poloidal_length`, `status`, `trace`, cumulative `L_par`/`L_pol`)
    - `wall`: list of wall contours used

    @param tokamaker TokaMaker object (current equilibrium is used) or TokaMaker equilibrium object
    @param sol_width Distance outside the separatrix at the OMP of the SOL surface used for connection lengths [m]
    @param wall Wall contour [n,2] or list of contours (default: limiter contours of the TokaMaker mesh)
    @param step Maximum poloidal step for field-line tracing [m]
    @param max_length Maximum poloidal length of a single trace [m] (default: 10x total wall length)
    @result Dictionary of divertor metrics (see above)
    '''
    if hasattr(tokamaker, '_tMaker_equil'):
        tMaker = tokamaker
        equil = tokamaker._tMaker_equil
        if equil is None:
            raise ValueError('TokaMaker object has no equilibrium')
    else:
        equil = tokamaker
        tMaker = tokamaker._tMaker
    if tMaker.settings.mirror_mode:
        raise ValueError('Divertor metrics are not supported in mirror mode')
    if wall is None:
        if tMaker.lim_contours is None:
            raise ValueError('No limiter contour available, specify "wall"')
        walls = tMaker.lim_contours
    elif isinstance(wall, numpy.ndarray) and wall.ndim == 2:
        walls = [wall]
    else:
        walls = list(wall)
    x_points, diverted = equil.get_xpoints()
    field = _TokaMakerField(equil)
    return _divertor_metrics(field, equil.o_point, x_points, walls, sol_width=sol_width, step=step,
                             max_length=max_length, diverted=bool(diverted))


def plot_divertor_metrics(fig, ax, metrics, wall_color='k', leg_color='tab:red', sol_color='tab:blue',
                          xpoint_color='k', annotate=True):
    r'''! Plot strike points, separatrix legs and the SOL field line from @ref OpenFUSIONToolkit.TokaMaker.divertor.compute_divertor_metrics "compute_divertor_metrics"

    @param fig Figure to add to
    @param ax Axis to add to
    @param metrics Output of @ref OpenFUSIONToolkit.TokaMaker.divertor.compute_divertor_metrics "compute_divertor_metrics"
    @param wall_color Color for wall contours (`None` to disable)
    @param leg_color Color for separatrix legs and strike points
    @param sol_color Color for the SOL field line
    @param xpoint_color Color for X-point markers (primary filled, secondary open)
    @param annotate Label each target with its name and \f$ f_x \f$
    '''
    if wall_color is not None:
        for contour in metrics['wall']:
            ax.plot(contour[:, 0], contour[:, 1], color=wall_color, lw=1.0)
    for end in metrics['sol']['ends']:
        ax.plot(end['trace'][:, 0], end['trace'][:, 1], color=sol_color, lw=1.0)
    for name, target in metrics['targets'].items():
        ax.plot(target['trace'][:, 0], target['trace'][:, 1], color=leg_color, lw=1.5,
                ls='-' if target['primary'] else '--')
        if target['status'] == 'wall':
            ax.plot(target['R_target'], target['Z_target'], 'o', color=leg_color)
            if annotate:
                ax.annotate('{0}\n$f_x$={1:.1f}'.format(name, target['poloidal_flux_expansion']),
                            (target['R_target'], target['Z_target']), fontsize=7,
                            textcoords='offset points', xytext=(4, 4))
    for xpt in metrics['xpoints']:
        if xpt['inside_wall']:
            ax.plot(xpt['R'], xpt['Z'], 'x' if xpt['primary'] else '+', color=xpoint_color, ms=8)
    ax.plot(metrics['R_OMP'], metrics['Z_OMP'], 's', color=sol_color, ms=4)
    ax.set_aspect('equal')
