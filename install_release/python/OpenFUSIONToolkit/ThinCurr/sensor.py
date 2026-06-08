#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! ThinCurr sensor definition and manipulation functionality

@authors Chris Hansen
@date May 2024
@ingroup doxy_oft_python
'''
import numpy


class flux_loop():
    '''! ThinCurr flux loop class'''
    def __init__(self, pts, name, scale=1.0):
        '''! Initialize flux loop object

        @param pts Points forming flux loop [:,3]
        @param name Name of flux loop
        @param scale Scale factor for signal
        '''
        self._pts = numpy.asarray(pts)
        if pts.shape[1] != 3:
            raise ValueError("Point array must be [:,3]")
        self._name = name
        self._scale = scale

    def write(self,fid):
        '''! Write flux loop to file

        @param fid File object for writing
        '''
        fid.write('\n{0} {1} {2}\n'.format(self._pts.shape[0], self._scale, self._name))
        for i in range(self._pts.shape[0]):
            fid.write('{0:.6E} {1:.6E} {2:.6E}\n'.format(*self._pts[i,:]))


class circular_flux_loop(flux_loop):
    '''! ThinCurr circular flux loop class'''
    def __init__(self, R, Z, name, scale=1.0, npts=180):
        '''! Initialize circular flux loop object

        @param R Radial position of flux loop
        @param Z Vertical position of flux loop
        @param name Name of flux loop
        @param scale Scale factor for signal
        @param npts Number of points used to define flux loop
        '''
        self._pts = numpy.zeros((npts,3))
        theta = numpy.linspace(0.0,2.0*numpy.pi,npts)
        self._pts[:,0] = R*numpy.cos(theta)
        self._pts[:,1] = R*numpy.sin(theta)
        self._pts[:,2] = Z
        self._name = name
        self._scale = scale


class Mirnov(flux_loop):
    '''! ThinCurr Mirnov sensor class'''
    def __init__(self, pt, norm, name, dx=1.E-5):
        '''! Initialize Mirnov sensor object

        @param pt Location of Mirnov sensor [3]
        @param norm Normal direction of Mirnov sensor [3]
        @param name Name of Mirnov sensor
        @param dx Scale of Mirnov flux loop
        '''
        center = numpy.asarray(pt)
        if center.shape[0] != 3:
            raise ValueError('"center" must be a 3-dimensional vector')
        norm = numpy.asarray(norm)
        if norm.shape[0] != 3:
            raise ValueError('"norm" must be a 3-dimensional vector')
        norm = norm/numpy.linalg.norm(norm)
        v1 = numpy.cross([1.0, 0.0, 0.0], norm)
        if numpy.linalg.norm(v1) < 1.E-8:
            v1 = numpy.cross([0.0, 1.0, 0.0], norm)
        v1 = v1/numpy.linalg.norm(v1)
        v2 = numpy.cross(v1, norm)
        v2 /= numpy.linalg.norm(v2)
        #
        self._pts = numpy.array([
            center+dx*(v1+v2)/2.0,
            center+dx*(-v1+v2)/2.0,
            center+dx*(-v1-v2)/2.0,
            center+dx*(v1-v2)/2.0,
            center+dx*(v1+v2)/2.0
        ])
        self._name = name
        self._scale = 1.0/numpy.power(dx, 2)


def save_sensors(sensor_list,filename='floops.loc'):
    '''! Save ThinCurr sensors to file

    @param sensor_list List of sensors
    @param filename Filename for output file
    '''
    for sensor in sensor_list:
        if not isinstance(sensor, flux_loop):
            raise ValueError("Invalid sensor object")
    with open(filename, 'w+') as fid:
        fid.write('{0}\n'.format(len(sensor_list)))
        for sensor in sensor_list:
            sensor.write(fid)