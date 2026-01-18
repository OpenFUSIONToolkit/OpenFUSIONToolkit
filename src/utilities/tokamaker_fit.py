# Constraint definition and handling for the Open FUSION Toolkit (OFT)
import numpy as np

field_con_id = 1
iTor_con_id = 2
eLam_con_id = 4
fluxLoop_con_id = 7
dFlux_con_id = 8
Pe_con_id = 9
q_con_id = 10
saddle_con_id = 11


class field_con:
    def __init__(self, pt=None, phi=0., norm=None, val=None, err=None):
        self.pt = pt
        self.phi = phi
        self.norm = norm
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        self.phi = float(values[2])
        values = file.readline().split()
        self.norm = (float(values[0]), float(values[1]), float(values[2]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(field_con_id))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.pt[0], self.pt[1], self.phi))
        file.write(' {0:E} {1:E} {2:E}\n'.format(self.norm[0], self.norm[1], self.norm[2]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class iTor_con:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(iTor_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class eLam_con:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(eLam_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class fluxLoop_con:
    def __init__(self, pt=None, val=None, err=None):
        self.pt = pt
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(fluxLoop_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt[0], self.pt[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class dFlux_con:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(dFlux_con_id))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class Pe_con:
    def __init__(self, pt=None, val=None, err=None):
        self.pt = pt
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(Pe_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt[0], self.pt[1]))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class q_con:
    def __init__(self, type=None, val=None, err=None, loc=0.):
        self.type = type
        self.val = val
        self.err = err
        self.loc = loc

    def read(self, file):
        values = file.readline().split()
        self.type = int(values[0])
        self.loc = float(values[1])
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(q_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.type, self.loc))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


class saddle_con:
    def __init__(self, pt1=None, pt2=None, width=None, val=None, err=None):
        self.pt1 = pt1
        self.pt2 = pt2
        self.width = width
        self.val = val
        self.err = err

    def read(self, file):
        values = file.readline().split()
        self.pt1 = (float(values[0]), float(values[1]))
        values = file.readline().split()
        self.pt2 = (float(values[0]), float(values[1]))
        value = file.readline()
        self.width = float(value)
        values = file.readline().split()
        self.val = float(values[0])
        self.err = 1./float(values[1])

    def write(self, file):
        file.write('{0}\n'.format(saddle_con_id))
        file.write(' {0:E} {1:E}\n'.format(self.pt1[0], self.pt1[1]))
        file.write(' {0:E} {1:E}\n'.format(self.pt2[0], self.pt2[1]))
        file.write(' {0:E}\n'.format(self.width))
        file.write(' {0:E} {1:E}\n\n'.format(self.val, 1./self.err))


con_map = {
    field_con_id: field_con,
    iTor_con_id: iTor_con,
    eLam_con_id: eLam_con,
    fluxLoop_con_id: fluxLoop_con,
    dFlux_con_id: dFlux_con,
    Pe_con_id: Pe_con,
    q_con_id: q_con,
    saddle_con_id: saddle_con
}


def write_fit_in(constraints, filename='fit.in'):
    ncons = len(constraints)
    with open(filename, 'w+') as fid:
        fid.write('{0:d}\n\n'.format(ncons))
        for con in constraints:
            con.write(fid)


def read_fit_in(filename='fit.in'):
    constraints = []
    with open(filename, 'r') as fid:
        ncons = int(fid.readline())
        for i in range(ncons):
            fid.readline()
            con_type = int(fid.readline())
            new_con_class = con_map[con_type]
            new_con = new_con_class()
            new_con.read(fid)
            constraints.append((i, new_con))
    return constraints


class coil_group:
    def __init__(self, center, size, nlayers, symmetric, current, name=None):
        self.center = center
        self.size = size
        self.nlayers = nlayers
        self.symmetric = symmetric
        self.current = current
        self.name = name

    def write(self, file, indent='  '):
        # Setup coil sizes
        dx = self.size[0]/self.nlayers[0]
        x0 = self.center[0] - dx*(self.nlayers[0]/2. - .5)
        dy = self.size[1]/self.nlayers[1]
        y0 = self.center[1] - dy*(self.nlayers[1]/2. - .5)
        # Create coils
        if self.name is None:
            string = indent + '<coil_set current="{0:E}">\n'.format(self.current)
        else:
            string = indent + '<coil_set current="{0:E}" name="{1}">\n'.format(self.current, self.name)
        for i in range(self.nlayers[0]):
            xc = x0 + i*dx
            for j in range(self.nlayers[1]):
                yc = y0 + j*dy
                string = string + indent + "  <coil>{0:f}, {1:f}</coil>\n".format(xc, yc)
                if self.symmetric:
                    string = string + indent + "  <coil>{0:f}, {1:f}</coil>\n".format(xc, -yc)
        # Close coil set group
        string = string + indent + '</coil_set>\n'
        file.write(string)


class coil_region:
    def __init__(self, id, current, vcont_gain=None, name=None):
        self.id = id
        self.current = current
        self.vcont_gain = vcont_gain
        self.name = name

    def write(self, file, indent='  '):
        # Create coils
        if self.name is None:
            string = indent + '<region type="coil" id="{0}">\n'.format(self.id)
        else:
            string = indent + '<region type="coil" id="{0}" name="{1}">\n'.format(self.id, self.name)
        string = string + indent + '  <current>{0}</current>\n'.format(self.current)
        if self.vcont_gain is not None:
            string = string + indent + '  <vcont_gain>{0}</vcont_gain>\n'.format(self.vcont_gain)
        string = string + indent + '</region>\n'
        file.write(string)


class cond_region:
    def __init__(self, id, neigs=0, weights=None, contiguous=True, limiter=True):
        self.id = id
        self.neigs = neigs
        self.weights = weights
        self.contiguous = contiguous
        self.limiter = limiter

    def write(self, file, indent='  '):
        # Build weight list
        weight_list = self.weights
        if weight_list is None:
            weight_list = np.zeros((self.neigs,))
        first = True
        weights = ''
        for weight in weight_list:
            if not first:
                weights = weights + ', '
            weights = weights + '{0:E}'.format(weight)
            first = False
        # Add XML entry
        string = indent + '<region type="wall" id="{0}">\n'.format(self.id)
        string = string + indent + '  <neigs>{0}</neigs>\n'.format(self.neigs)
        string = string + indent + '  <weights>{0}</weights>\n'.format(weights)
        if not self.contiguous:
            string = string + indent + '  <contiguous>false</contiguous>\n'
        if not self.limiter:
            string = string + indent + '  <limiter>false</limiter>\n'
        string = string + indent + '</region>\n'
        file.write(string)


def read_eqdsk(filename):
    def read_1d(fid, j, n):
        output = np.zeros((n,))
        for i in range(n):
            if j == 0:
                line = fid.readline()
            output[i] = line[j:j+16]
            j += 16
            if j == 16*5:
                j = 0
        return output, j

    def read_2d(fid, j, n, m):
        output = np.zeros((m, n))
        for k in range(n):
            for i in range(m):
                if j == 0:
                    line = fid.readline()
                output[i, k] = line[j:j+16]
                j += 16
                if j == 16*5:
                    j = 0
        return output, j
    # Read-in data
    eqdsk_obj = {}
    with open(filename, 'r') as fid:
        # Get sizes
        line = fid.readline()
        eqdsk_obj['case'] = line[:48]
        split_line = line[48:].split()
        eqdsk_obj['nz'] = int(split_line[-1])
        eqdsk_obj['nr'] = int(split_line[-2])
        # Read header content
        line_keys = [['rdim',  'zdim',  'raxis',  'rleft',  'zmid'],
                     ['raxis', 'zaxis', 'psimax', 'psimin', 'bcentr'],
                     ['itor',  'skip',  'skip',   'skip',   'skip'],
                     ['skip',  'skip',  'skip',   'skip',   'skip']]
        for i in range(4):
            line = fid.readline()
            for j in range(5):
                if line_keys[i][j] == 'skip':
                    continue
                line_seg = line[j*16:(j+1)*16]
                eqdsk_obj[line_keys[i][j]] = float(line_seg)
        # Read flux profiles
        j = 0
        keys = ['fpol', 'pres', 'ffprim', 'pprime']
        for key in keys:
            eqdsk_obj[key], j = read_1d(fid, j, eqdsk_obj['nr'])
        # Read PSI grid
        eqdsk_obj['psirz'], j = read_2d(fid, j, eqdsk_obj['nz'],
                                        eqdsk_obj['nr'])
        # Read q-profile
        eqdsk_obj['qpsi'], j = read_1d(fid, j, eqdsk_obj['nr'])
        # Read limiter count
        line = fid.readline()
        eqdsk_obj['nlim'] = int(line.split()[1])
        # Read outer flux surface
        eqdsk_obj['rzout'], j = read_2d(fid, j, eqdsk_obj['nr'], 2)
        # Read limiting corners
        eqdsk_obj['rzlim'], j = read_2d(fid, j, eqdsk_obj['nlim'], 2)
    return eqdsk_obj

def write_eqdsk(eqdsk_in, filename):
    def write_1d(fid, j, input):
        n = input.shape[0]
        for i in range(n):
            if j == 0:
                fid.write('\n')
            fid.write('{0:16.9E}'.format(input[i]))
            j += 16
            if j == 16*5:
                j = 0
        return j

    def write_2d(fid, j, input):
        m = input.shape[0]
        n = input.shape[1]
        for k in range(n):
            for i in range(m):
                if j == 0:
                    fid.write('\n')
                fid.write('{0:16.9E}'.format(input[i, k]))
                j += 16
                if j == 16*5:
                    j = 0
        return j
    # Read-in data
    eqdsk_obj = eqdsk_in.copy()
    eqdsk_obj['xdum'] = 0.0
    with open(filename, 'w+') as fid:
        # Write sizes
        fid.write(eqdsk_obj['case'])
        fid.write('{0:4d}{1:4d}{2:4d}'.format(0, eqdsk_obj['nz'], eqdsk_obj['nr']))
        # Write header content
        line_keys = [['rdim',  'zdim',   'raxis',  'rleft',  'zmid'],
                     ['raxis', 'zaxis',  'psimax', 'psimin', 'bcentr'],
                     ['itor',  'psimax', 'xdum',   'raxis',  'xdum'],
                     ['zaxis', 'xdum',   'psimin', 'xdum',   'xdum']]
        for i in range(4):
            fid.write('\n')
            for j in range(5):
                fid.write('{0:16.9E}'.format(eqdsk_obj[line_keys[i][j]]))
        # Write flux profiles
        j = 0
        keys = ['fpol', 'pres', 'ffprim', 'pprime']
        for key in keys:
            j = write_1d(fid, j, eqdsk_obj[key])
        # Write PSI grid
        j = write_2d(fid, j, eqdsk_obj['psirz'])
        # Write q-profile
        j = write_1d(fid, j, eqdsk_obj['qpsi'])
        # Write LCFS and limiter counts
        fid.write('\n{0:4d}{1:4d}'.format(eqdsk_obj['nr'], eqdsk_obj['nlim']))
        # Write outer flux surface
        j = write_2d(fid, j, eqdsk_obj['rzout'])
        # Write limiting corners
        j = write_2d(fid, j, eqdsk_obj['rzlim'])
    return eqdsk_obj
