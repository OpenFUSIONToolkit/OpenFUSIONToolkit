#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#------------------------------------------------------------------------------
#
# Python construction of H^1, H(Curl) elements for simplectic elements in 2D and
# 3D. Designed to be used for FE basis definition in the Open FUSION Toolkit (OFT)
#
#------------------------------------------------------------------------------
from __future__ import print_function
try:
    from sympy import *
except:
    print("================================")
    print("ERROR: SYMPY is required to run this script")
    print("================================")
    raise
else:
    from sympy.polys.orthopolys import legendre_poly as legendre
# Setup common symbolic variables
x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
s, t, x, zed = symbols('s t x zed')
#----------------------------------------------------------------
# Scaled Legendre Polynomial
#----------------------------------------------------------------
def ls(s,t,i):
    return (t**i)*legendre(i,s/t)
#----------------------------------------------------------------
# Scaled Integrated Legendre Polynomial
#----------------------------------------------------------------
def Ls(s,t,i):
    return (t**i)*integrate(legendre(i-1,x),(x,-1,s/t))
#----------------------------------------------------------------
# Basis functions for the H^1 space
#----------------------------------------------------------------
def H1(n):
    if n>0:
        fp, fe, ff, fc = H1(n-1)
    else:
        fp = []
        fe = []
        ff = []
        fc = []
        return fp, fe, ff, fc
    # Single point DOF
    if(n==1):
        fp.append(x1)
    # Print Edge DOFs
    for i in range(2,n+1):
        ui = simplify(Ls(x1-x2,x1+x2,i))
        f = simplify(ui)
        if f in fe:
            continue
        fe.append(f)
    # Print Face DOFs
    for i in range(2,n):
        for j in range(n-i):
            ui = simplify(Ls(x1-x2,x1+x2,i))
            vj = simplify(x3*ls(x3-x1-x2,x1+x2+x3,j))
            f = simplify(vj*ui)
            if f in ff:
                continue
            ff.append(f)
    # Print Cell DOFs
    for i in range(2,n-1):
        for j in range(n-i-1):
            for k in range(n-1-i-j):
                ui = simplify(Ls(x1-x2,x1+x2,i))
                vj = simplify(x3*ls(x3-x1-x2,x1+x2+x3,j))
                wk = simplify(x4*legendre(k,x4-x1-x2-x3))
                f = simplify(vj*wk*ui)
                if f in fc:
                    continue
                fc.append(f)
    return fp, fe, ff, fc
#----------------------------------------------------------------
# Basis functions for the H(Curl) space
#----------------------------------------------------------------
def HCurl(n):
    if n>1:
        fe, ff, fc = HCurl(n-1)
    else:
        fe = []
        ff = []
        fc = []
        return fe, ff, fc
    # Single edge DOF
    if(n==2):
        fe.append((x2,-x1))
    # Print Face DOFs
    for i in range(2,n):
        for j in range(n-i):
            ui = simplify(Ls(x1-x2,x1+x2,i))
            vj = simplify(x3*ls(x3-x1-x2,x1+x2+x3,j))
            f = (vj*ui.diff(x1) - ui*vj.diff(x1),vj*ui.diff(x2) - ui*vj.diff(x2),vj*ui.diff(x3) - ui*vj.diff(x3))
            if not(f in ff):
                ff.append(f)
            if (j<n):
                f = (x2*vj,-x1*vj,zed)
                if f in ff:
                    continue
                ff.append(f)
    # Print Cell DOFs
    for i in range(2,n-1):
        for j in range(n-i-1):
            for k in range(n-1-i-j):
                ui = simplify(Ls(x1-x2,x1+x2,i))
                vj = simplify(x3*ls(x3-x1-x2,x1+x2+x3,j))
                wk = simplify(x4*legendre(k,x4-x1-x2-x3))
                f=(vj*wk*ui.diff(x1) - ui*wk*vj.diff(x1) + ui*vj*wk.diff(x1),vj*wk*ui.diff(x2) - ui*wk*vj.diff(x2) + ui*vj*wk.diff(x2),vj*wk*ui.diff(x3) - ui*wk*vj.diff(x3) + ui*vj*wk.diff(x3),vj*wk*ui.diff(x4) - ui*wk*vj.diff(x4) + ui*vj*wk.diff(x4))
                if not(f in fc):
                    fc.append(f)
                    f=(vj*wk*ui.diff(x1) - ui*wk*vj.diff(x1) - ui*vj*wk.diff(x1),vj*wk*ui.diff(x2) - ui*wk*vj.diff(x2) - ui*vj*wk.diff(x2),vj*wk*ui.diff(x3) - ui*wk*vj.diff(x3) - ui*vj*wk.diff(x3),vj*wk*ui.diff(x4) - ui*wk*vj.diff(x4) - ui*vj*wk.diff(x4))
                    fc.append(f)
                if (j+k<n-1):
                    f=(x2*vj*wk,-x1*vj*wk,zed,zed)
                    if f in fc:
                        continue
                    fc.append(f)
    return fe, ff, fc
#----------------------------------------------------------------
# Basis functions for the H(Div) space, ie. Curl( H(Curl) )
#----------------------------------------------------------------
def HCurlCurl(n):
    # Return if invalid order
    if(n<2):
        return
    # Single edge DOF
    print("Edge DOF")
    print("g1xg2 = {0}".format(-2))
    # Print Face DOFs
    for i in range(2,n):
        for j in range(n-i):
            print("Face DOF Order = {0} {1}".format(i,j))
            ui = Ls(f1-f2,f1+f2,i)
            vj = f3*ls(f3-f1-f2,f1+f2+f3,j)
            print("Type 2")
            g1 = vj*ui.diff(f1) - ui*vj.diff(f1)
            g2 = vj*ui.diff(f2) - ui*vj.diff(f2)
            g3 = vj*ui.diff(f3) - ui*vj.diff(f3)
            print("g1xg2 = {0}".format(simplify(-g1.diff(f2) + g2.diff(f1))))
            print("g1xg3 = {0}".format(simplify(-g1.diff(f3) + g3.diff(f1))))
            print("g2xg3 = {0}".format(simplify(-g2.diff(f3) + g3.diff(f2))))
            if (j<n):
                print("Type 3")
                g1 = simplify(f2*vj)
                g2 = -simplify(f1*vj)
                g3 = 0*f1
                print("g1xg2 = {0}".format(simplify(-g1.diff(f2) + g2.diff(f1))))
                print("g1xg3 = {0}".format(simplify(-g1.diff(f3) + g3.diff(f1))))
                print("g2xg3 = {0}".format(simplify(-g2.diff(f3) + g3.diff(f2))))
    # Print Cell DOFs
    for i in range(2,n-1):
        for j in range(n-i-1):
            for k in range(n-1-i-j):
                print("Interior DOF Order = {0} {1} {2}".format(i,j,k))
                ui = Ls(f1-f2,f1+f2,i)
                vj = f3*ls(f3-f1-f2,f1+f2+f3,j)
                wk = f4*legendre(k,f4-f1-f2-f3)
                print("Type 2-1")
                g1 = vj*wk*ui.diff(f1) - ui*wk*vj.diff(f1) + ui*vj*wk.diff(f1)
                g2 = vj*wk*ui.diff(f2) - ui*wk*vj.diff(f2) + ui*vj*wk.diff(f2)
                g3 = vj*wk*ui.diff(f3) - ui*wk*vj.diff(f3) + ui*vj*wk.diff(f3)
                g4 = vj*wk*ui.diff(f4) - ui*wk*vj.diff(f4) + ui*vj*wk.diff(f4)
                print("g1xg2 = ", simplify(-g1.diff(f2) + g2.diff(f1)))
                print("g1xg3 = ", simplify(-g1.diff(f3) + g3.diff(f1)))
                print("g1xg4 = ", simplify(-g1.diff(f4) + g4.diff(f1)))
                print("g2xg3 = ", simplify(-g2.diff(f3) + g3.diff(f2)))
                print("g2xg4 = ", simplify(-g2.diff(f4) + g4.diff(f2)))
                print("g3xg4 = ", simplify(-g3.diff(f4) + g4.diff(f3)))
                print("Type 2-2")
                g1 = vj*wk*ui.diff(f1) - ui*wk*vj.diff(f1) - ui*vj*wk.diff(f1)
                g2 = vj*wk*ui.diff(f2) - ui*wk*vj.diff(f2) - ui*vj*wk.diff(f2)
                g3 = vj*wk*ui.diff(f3) - ui*wk*vj.diff(f3) - ui*vj*wk.diff(f3)
                g4 = vj*wk*ui.diff(f4) - ui*wk*vj.diff(f4) - ui*vj*wk.diff(f4)
                print("g1xg2 = ", simplify(-g1.diff(f2) + g2.diff(f1)))
                print("g1xg3 = ", simplify(-g1.diff(f3) + g3.diff(f1)))
                print("g1xg4 = ", simplify(-g1.diff(f4) + g4.diff(f1)))
                print("g2xg3 = ", simplify(-g2.diff(f3) + g3.diff(f2)))
                print("g2xg4 = ", simplify(-g2.diff(f4) + g4.diff(f2)))
                print("g3xg4 = ", simplify(-g3.diff(f4) + g4.diff(f3)))
                if (j+k<n-1):
                    print("Type 3")
                    g1 = f2*vj*wk
                    g2 = -f1*vj*wk
                    g3 = 0*f1
                    g4 = 0*f1
                    print("g1xg2 = ", simplify(-g1.diff(f2) + g2.diff(f1)))
                    print("g1xg3 = ", simplify(-g1.diff(f3) + g3.diff(f1)))
                    print("g1xg4 = ", simplify(-g1.diff(f4) + g4.diff(f1)))
                    print("g2xg3 = ", simplify(-g2.diff(f3) + g3.diff(f2)))
                    print("g2xg4 = ", simplify(-g2.diff(f4) + g4.diff(f2)))
                    print("g3xg4 = ", simplify(-g3.diff(f4) + g4.diff(f3)))
#----------------------------------------------------------------
# Script driver
#----------------------------------------------------------------
if __name__ == '__main__':
    #
    x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
    order = 1
    #
    import argparse
    parser = argparse.ArgumentParser(description='Generate conforming element basis functions for Tetrahedra.')
    parser.add_argument('-o', '--order', help='Desired interpolation order.', default=1)
    args=parser.parse_args()
    order=int(args.order)
    #
    print("\n====== Scalar Functions ======")
    fp, fe, ff, fc = H1(order)
    print("Point Functions")
    for i in range(len(fp)):
        f = fp[i]
        print(simplify(f))
    print("Edge Functions")
    for i in range(len(fe)):
        f = fe[i]
        print(i)
        print(simplify(f))
    print("Face Functions")
    for i in range(len(ff)):
        f = ff[i]
        print(i)
        print(simplify(f))
    print("Cell Functions")
    for i in range(len(fc)):
        f = fc[i]
        print(i)
        print(simplify(f))
    #
    print("\n====== Gradient Functions ======")
    print("Point Functions")
    for i in range(len(fp)):
        f = fp[i]
        print(i)
        print('g = {0}'.format(simplify(f.diff(x1))))
    print("Edge Functions")
    for i in range(len(fe)):
        f = fe[i]
        print(i)
        print('g1 = {0}'.format(simplify(f.diff(x1))))
        print('g2 = {0}'.format(simplify(f.diff(x2))))
    print("Face Functions")
    for i in range(len(ff)):
        f = ff[i]
        print(i)
        print('g1 = {0}'.format(simplify(f.diff(x1))))
        print('g2 = {0}'.format(simplify(f.diff(x2))))
        print('g3 = {0}'.format(simplify(f.diff(x3))))
    print("Cell Functions")
    for i in range(len(fc)):
        f = fc[i]
        print(i)
        print('g1 = {0}'.format(simplify(f.diff(x1))))
        print('g2 = {0}'.format(simplify(f.diff(x2))))
        print('g3 = {0}'.format(simplify(f.diff(x3))))
        print('g4 = {0}'.format(simplify(f.diff(x4))))
    #
    print("\n====== Second Derivatives ======")
    print("Point Functions")
    for i in range(len(fp)):
        f = fp[i]
        print(i)
        print('d =  {0}'.format(simplify(f.diff(x1).diff(x1))))
    print("Edge Functions")
    for i in range(len(fe)):
        f = fe[i]
        print(i)
        print('d1 = {0}'.format(simplify(f.diff(x1).diff(x1))))
        print('d2 = {0}'.format(simplify(f.diff(x1).diff(x2))))
        print('d3 = {0}'.format(simplify(f.diff(x2).diff(x2))))
    print("Face Functions")
    for i in range(len(ff)):
        f = ff[i]
        print(i)
        print('d1 = {0}'.format(simplify(f.diff(x1).diff(x1))))
        print('d2 = {0}'.format(simplify(f.diff(x1).diff(x2))))
        print('d3 = {0}'.format(simplify(f.diff(x1).diff(x3))))
        print('d4 = {0}'.format(simplify(f.diff(x2).diff(x2))))
        print('d5 = {0}'.format(simplify(f.diff(x2).diff(x3))))
        print('d6 = {0}'.format(simplify(f.diff(x3).diff(x3))))
    print("Cell Functions")
    for i in range(len(fc)):
        f = fc[i]
        print(i)
        print('d1 =  {0}'.format(simplify(f.diff(x1).diff(x1))))
        print('d2 =  {0}'.format(simplify(f.diff(x1).diff(x2))))
        print('d3 =  {0}'.format(simplify(f.diff(x1).diff(x3))))
        print('d4 =  {0}'.format(simplify(f.diff(x1).diff(x4))))
        print('d5 =  {0}'.format(simplify(f.diff(x2).diff(x2))))
        print('d6 =  {0}'.format(simplify(f.diff(x2).diff(x3))))
        print('d7 =  {0}'.format(simplify(f.diff(x2).diff(x4))))
        print('d8 =  {0}'.format(simplify(f.diff(x3).diff(x3))))
        print('d9 =  {0}'.format(simplify(f.diff(x3).diff(x4))))
        print('d10 = {0}'.format(simplify(f.diff(x4).diff(x4))))
    #
    print("\n====== Curl Functions ======")
    fe, ff, fc = HCurl(order)
    print("Edge Functions")
    for i in range(len(fe)):
        f = fe[i]
        print(i)
        print('g1 = {0}'.format(simplify(f[0].evalf(subs={zed:0}))))
        print('g2 = {0}'.format(simplify(f[1].evalf(subs={zed:0}))))
    print("Face Functions")
    for i in range(len(ff)):
        f = ff[i]
        print(i)
        print('g1 = {0}'.format(simplify(f[0].evalf(subs={zed:0}))))
        print('g2 = {0}'.format(simplify(f[1].evalf(subs={zed:0}))))
        print('g3 = {0}'.format(simplify(f[2].evalf(subs={zed:0}))))
    print("Cell Functions")
    for i in range(len(fc)):
        f = fc[i]
        print(i)
        print('g1 = {0}'.format(simplify(f[0].evalf(subs={zed:0}))))
        print('g2 = {0}'.format(simplify(f[1].evalf(subs={zed:0}))))
        print('g3 = {0}'.format(simplify(f[2].evalf(subs={zed:0}))))
        print('g4 = {0}'.format(simplify(f[3].evalf(subs={zed:0}))))
    #
    print("\n====== Curl(Curl) Functions ======")
    print("Edge Functions")
    for i in range(len(fe)):
        f = fe[i]
        print(i)
        print("g1xg2 = {0}".format(simplify((-f[0].diff(x2) + f[1].diff(x1)).evalf(subs={zed:0}))))
    print("Face Functions")
    for i in range(len(ff)):
        f = ff[i]
        print(i)
        print("g1xg2 = {0}".format(simplify((-f[0].diff(x2) + f[1].diff(x1)).evalf(subs={zed:0}))))
        print("g1xg3 = {0}".format(simplify((-f[0].diff(x3) + f[2].diff(x1)).evalf(subs={zed:0}))))
        print("g2xg3 = {0}".format(simplify((-f[1].diff(x3) + f[2].diff(x2)).evalf(subs={zed:0}))))
    print("Cell Functions")
    for i in range(len(fc)):
        f = fc[i]
        print(i)
        print("g1xg2 = {0}".format(simplify((-f[0].diff(x2) + f[1].diff(x1)).evalf(subs={zed:0}))))
        print("g1xg3 = {0}".format(simplify((-f[0].diff(x3) + f[2].diff(x1)).evalf(subs={zed:0}))))
        print("g1xg4 = {0}".format(simplify((-f[0].diff(x4) + f[3].diff(x1)).evalf(subs={zed:0}))))
        print("g2xg3 = {0}".format(simplify((-f[1].diff(x3) + f[2].diff(x2)).evalf(subs={zed:0}))))
        print("g2xg4 = {0}".format(simplify((-f[1].diff(x4) + f[3].diff(x2)).evalf(subs={zed:0}))))
        print("g3xg4 = {0}".format(simplify((-f[2].diff(x4) + f[3].diff(x3)).evalf(subs={zed:0}))))
