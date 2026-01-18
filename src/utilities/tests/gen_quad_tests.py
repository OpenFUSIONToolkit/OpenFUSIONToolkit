#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#------------------------------------------------------------------------------
#
# Class definition for manipulating Open FUSION Toolkit (OFT) structured binary files in Python
#
#------------------------------------------------------------------------------
from __future__ import print_function
try:
    import numpy
except:
    print("================================")
    print("ERROR: NUMPY is required to run this script")
    print("================================")
    raise
try:
    from sympy import symbols, Function, Integral
except:
    print("================================")
    print("ERROR: SYMPY is required to run this script")
    print("================================")
    raise
# Setup common symbolic variables
x, y, z = symbols('x y z')
f = symbols('f', cls=Function)
#----------------------------------------------------------------
# Get exponents for all combinations up to a given order in 1D
#----------------------------------------------------------------
def get_exponents_1d(order):
    e = numpy.array([[0,]])
    for i in range(1,order):
        e = numpy.append([[i,]],e,axis=0)
    return e
#----------------------------------------------------------------
# Get exponents for all combinations up to a given order in 2D
#----------------------------------------------------------------
def get_exponents_2d(order):
    e = numpy.array([[1,0]])
    for i in range(order):
        for j in range(order):
            if i+j>order:
                continue
            e = numpy.append([[i,j]],e,axis=0)
    return e
#----------------------------------------------------------------
# Get exponents for all combinations up to a given order in 3D
#----------------------------------------------------------------
def get_exponents_3d(order):
    e = numpy.array([[1,0,0]])
    for i in range(order):
        for j in range(order):
            for k in range(order):
                if i+j+k>order:
                    continue
                e = numpy.append([[i,j,k]],e,axis=0)
    return e
#----------------------------------------------------------------
# Script driver
#----------------------------------------------------------------
if __name__ == '__main__':
    # Create quadrature tests for 1D rules
    e = get_exponents_1d(8)
    fid = open('quad_1d.tests','w+')
    fid.write('{0:12d}\n'.format(e.shape[0]))
    for i in range(e.shape[0]):
        f = x**e[i,0]
        expected = Integral(f,(x,0,1)).doit()
        es = (e[i,0],)
        fid.write('{0:12d}{1:25.17E}{2:25.17E}\n'.format(es[0],1,float(expected)))
    fid.close()
    # Create quadrature tests for 2D rules
    e = get_exponents_2d(19)
    fid = open('quad_2d.tests','w+')
    fid.write('{0:12d}\n'.format(e.shape[0]))
    for i in range(e.shape[0]):
        f = x**e[i,0] + y**e[i,1]
        expected = Integral(f,(x,0,1),(y,0,1)).doit()
        es = (e[i,0],e[i,1])
        fid.write('{0:12d}{1:12d}{2:25.17E}{3:25.17E}{4:25.17E}\n'.format(es[0],es[1],1,1,float(expected)))
    fid.close()
    # Create quadrature tests for 3D rules
    e = get_exponents_3d(12)
    fid = open('quad_3d.tests','w+')
    fid.write('{0:12d}\n'.format(e.shape[0]))
    for i in range(e.shape[0]):
        f = x**e[i,0] + y**e[i,1] + z**e[i,2]
        expected = Integral(f,(x,0,1),(y,0,1),(z,0,1)).doit()
        es = (e[i,0],e[i,1],e[i,2])
        fid.write('{0:12d}{1:12d}{2:12d}{3:25.17E}{4:25.17E}{5:25.17E}{6:25.17E}\n'.format(es[0],es[1],es[2],1,1,1,float(expected)))
    fid.close()
