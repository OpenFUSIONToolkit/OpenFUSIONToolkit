#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#------------------------------------------------------------------------------
#
# Python construction of Lagrange elements for simplectic elements in 2D and
# 3D. Designed to be used for FE basis definition in the Open FUSION Toolkit (OFT)
#
#------------------------------------------------------------------------------
from __future__ import print_function
import sys
try:
    import numpy
except:
    print("================================")
    print("ERROR: NUMPY is required to run this script")
    print("================================")
    raise
try:
    from sympy import *
except:
    print("================================")
    print("ERROR: SYMPY is required to run this script")
    print("================================")
    raise
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_nodes(order):
        eds = []
        for i in range(order-1):
            eds.append((i+1)*1./(order))
        return eds
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_permute2(order):
    sets = []
    for i in range(order+1):
        j = order-i
        sets.append((i,j))
    return sets
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_permute3(order):
    sets = []
    for i in range(order+1):
        for j in range(order+1):
            if i+j > order:
                break
            k = order-i-j
            sets.append((i,j,k))
    return sets
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_permute4(order):
    sets = []
    for i in range(order+1):
        for j in range(order+1):
            if i+j > order:
                break
            for k in range(order+1):
                if i+j+k > order:
                    break
                l = order-i-j-k
                sets.append((i,j,k,l))
    return sets
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_fp(order):
    x1 = Symbol('x1')
    f = Function('f')
    nodes = get_nodes(order)
    f = x1
    node_loc = 1.
    for node in nodes:
        f = f*(x1-node)
    norm = f.evalf(subs={x1:1})
    return f/norm, node_loc
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_fe(order):
    if order < 2:
        return None
    x1, x2 = symbols('x1 x2')
    f = Function('f')
    funs = []
    node_locs = []
    sets = get_permute2(order-2)
    nodes = get_nodes(order)
    for myset in sets:
        f = x1*x2
        for j in range(myset[0]):
            f = f*(x1-nodes[j])
        for j in range(myset[1]):
            f = f*(x2-nodes[j])
        norm = f.evalf(subs={x1:nodes[myset[0]], x2:nodes[myset[1]]})
        funs.append(f/norm)
        node_locs.append((nodes[myset[0]],nodes[myset[1]]))
    return funs, node_locs
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_ff(order):
    if order < 3:
        return None
    x1, x2, x3 = symbols('x1 x2 x3')
    f = Function('f')
    funs = []
    node_locs = []
    sets = get_permute3(order-3)
    nodes = get_nodes(order)
    for myset in sets:
        f = x1*x2*x3
        for j in range(myset[0]):
            f = f*(x1-nodes[j])
        for j in range(myset[1]):
            f = f*(x2-nodes[j])
        for j in range(myset[2]):
            f = f*(x3-nodes[j])
        norm = f.evalf(subs={x1:nodes[myset[0]], x2:nodes[myset[1]], x3:nodes[myset[2]]})
        funs.append(f/norm)
        node_locs.append((nodes[myset[0]],nodes[myset[1]],nodes[myset[2]]))
    return funs, node_locs
#----------------------------------------------------------------
#
#----------------------------------------------------------------
def get_fc(order):
    if order < 4:
        return None
    x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
    f = Function('f')
    funs = []
    node_locs = []
    sets = get_permute4(order-4)
    nodes = get_nodes(order)
    for myset in sets:
        f = x1*x2*x3*x4
        for j in range(myset[0]):
            f = f*(x1-nodes[j])
        for j in range(myset[1]):
            f = f*(x2-nodes[j])
        for j in range(myset[2]):
            f = f*(x3-nodes[j])
        for j in range(myset[3]):
            f = f*(x4-nodes[j])
        norm = f.evalf(subs={x1:nodes[myset[0]], x2:nodes[myset[1]], x3:nodes[myset[2]], x4:nodes[myset[3]]})
        funs.append(f/norm)
        node_locs.append((nodes[myset[0]],nodes[myset[1]],nodes[myset[2]],nodes[myset[3]]))
    return funs, node_locs
#----------------------------------------------------------------
#
#----------------------------------------------------------------
class lagrange_interp:
    def __init__(self,order=1):
        self.funs_point, self.nodes_point = get_fp(order)
        self.np = 1
        if order>1:
            self.funs_edge, self.nodes_edge = get_fe(order)
            self.ne = len(self.funs_edge)
        else:
            self.funs_edge = []
            self.nodes_edge = []
            self.ne = 0
        if order>2:
            self.funs_face, self.nodes_face = get_ff(order)
            self.nf = len(self.funs_face)
        else:
            self.funs_face = []
            self.nodes_face = []
            self.nf = 0
        if order>3:
            self.funs_cell, self.nodes_cell = get_fc(order)
            self.nc = len(self.funs_cell)
        else:
            self.funs_cell = []
            self.nodes_cell = []
            self.nc = 0

    def eval_point(self,u):
        x1 = symbols('x1')
        return self.funs_point.evalf(subs={x1:u})

    def eval_edge(self,u,i):
        x1, x2 = symbols('x1 x2')
        return self.funs_edge[i].evalf(subs={x1:u[0], x2:u[1]})

    def eval_face(self,u,i):
        x1, x2, x3 = symbols('x1 x2 x3')
        return self.funs_face[i].evalf(subs={x1:u[0], x2:u[1], x3:u[2]})

    def eval_cell(self,u,i):
        x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
        return self.funs_cell[i].evalf(subs={x1:u[0], x2:u[1], x3:u[2], x4:u[3]})
#----------------------------------------------------------------
# Script driver
#----------------------------------------------------------------
if __name__ == '__main__':
    #
    x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
    order = 1
    #
    import argparse
    parser = argparse.ArgumentParser(description='Generate Lagrange basis functions for Tetrahedra.')
    parser.add_argument('-o', '--order', help='Desired interpolation order.', default=1)
    args=parser.parse_args()
    order=int(args.order)
    #
    print('Point DOF')
    funs_point, nodes_point = get_fp(order)
    print('f = ', funs_point)
    #
    if order > 1:
        print('Edge DOF')
        funs_edge, nodes_edge = get_fe(order)
        for i in range(len(funs_edge)):
            print(i)
            print('f = ', funs_edge[i])
    #
    if order > 2:
        print('Face DOF')
        funs_face, nodes_face = get_ff(order)
        for i in range(len(funs_face)):
            print(i)
            print('f = ', funs_face[i])
    #
    if order > 3:
        print('Cell DOF')
        funs_cell, nodes_cell = get_fc(order)
        for i in range(len(funs_cell)):
            print(i)
            print('f = ', funs_cell[i])

    #
    print('Point DOF')
    print('g = ', simplify(funs_point.diff(x1)))
    #
    if order > 1:
        print('Edge DOF')
        for i in range(len(funs_edge)):
            print(i)
            print('g1 = ', simplify(funs_edge[i].diff(x1)))
            print('g2 = ', simplify(funs_edge[i].diff(x2)))
    #
    if order > 2:
        print('Face DOF')
        for i in range(len(funs_face)):
            print(i)
            print('g1 = ', simplify(funs_face[i].diff(x1)))
            print('g2 = ', simplify(funs_face[i].diff(x2)))
            print('g3 = ', simplify(funs_face[i].diff(x3)))
    #
    if order > 3:
        print('Cell DOF')
        for i in range(len(funs_cell)):
            print(i)
            print('g1 = ', simplify(funs_cell[i].diff(x1)))
            print('g2 = ', simplify(funs_cell[i].diff(x2)))
            print('g3 = ', simplify(funs_cell[i].diff(x3)))
            print('g4 = ', simplify(funs_cell[i].diff(x4)))
    #
    print('Point DOF')
    print('d = ', simplify(funs_point.diff(x1).diff(x1)))
    #
    if order > 1:
        print('Edge DOF')
        for i in range(len(funs_edge)):
            print(i)
            print('d1 = ', simplify(funs_edge[i].diff(x1).diff(x1)))
            print('d2 = ', simplify(funs_edge[i].diff(x1).diff(x2)))
            print('d3 = ', simplify(funs_edge[i].diff(x2).diff(x2)))
    #
    if order > 2:
        print('Face DOF')
        for i in range(len(funs_face)):
            print(i)
            print('d1 = ', simplify(funs_face[i].diff(x1).diff(x1)))
            print('d2 = ', simplify(funs_face[i].diff(x1).diff(x2)))
            print('d3 = ', simplify(funs_face[i].diff(x1).diff(x3)))
            print('d4 = ', simplify(funs_face[i].diff(x2).diff(x2)))
            print('d5 = ', simplify(funs_face[i].diff(x2).diff(x3)))
            print('d6 = ', simplify(funs_face[i].diff(x3).diff(x3)))
    #
    if order > 3:
        print('Cell DOF')
        for i in range(len(funs_cell)):
            print(i)
            print('d1 = ', simplify(funs_cell[i].diff(x1).diff(x1)))
            print('d2 = ', simplify(funs_cell[i].diff(x1).diff(x2)))
            print('d3 = ', simplify(funs_cell[i].diff(x1).diff(x3)))
            print('d4 = ', simplify(funs_cell[i].diff(x1).diff(x4)))
            print('d5 = ', simplify(funs_cell[i].diff(x2).diff(x2)))
            print('d6 = ', simplify(funs_cell[i].diff(x2).diff(x3)))
            print('d7 = ', simplify(funs_cell[i].diff(x2).diff(x4)))
            print('d8 = ', simplify(funs_cell[i].diff(x3).diff(x3)))
            print('d9 = ', simplify(funs_cell[i].diff(x3).diff(x4)))
            print('d10 = ', simplify(funs_cell[i].diff(x4).diff(x4)))

    #
    print('Point Nodes')
    print('n = ', nodes_point)
    #
    if order > 1:
        print('Edge Nodes')
        for i in range(len(nodes_edge)):
            print(i)
            print('n = ', nodes_edge[i])
    #
    if order > 2:
        print('Face Nodes')
        for i in range(len(nodes_face)):
            print(i)
            print('n = ', nodes_face[i])
    #
    if order > 3:
        print('Cell Nodes')
        for i in range(len(nodes_cell)):
            print(i)
            print('n = ', nodes_cell[i])
