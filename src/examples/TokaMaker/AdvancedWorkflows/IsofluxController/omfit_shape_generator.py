#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 09:56:11 2025

@author: mparsons

These functions are taken directly from OMFIT

"""

import numpy as np
import contourpy

##################################################################

def boundaryShape(
    a,
    eps,
    kapu,
    kapl,
    delu,
    dell,
    zetaou,
    zetaiu,
    zetail,
    zetaol,
    zoffset,
    upnull=False,
    lonull=False,
    npts=90,
    doPlot=False,
    newsq=np.zeros(4),
    **kw,
):
    '''
    Function used to generate boundary shapes based on `T. C. Luce, PPCF, 55 9 (2013)`
    Direct Python translation of the IDL program /u/luce/idl/shapemaker3.pro

    :param a: minor radius

    :param eps: aspect ratio

    :param kapu: upper elongation

    :param lkap: lower elongation

    :param delu: upper triangularity

    :param dell: lower triangularity

    :param zetaou: upper outer squareness

    :param zetaiu: upper inner squareness

    :param zetail: lower inner squareness

    :param zetaol: lower outer squareness

    :param zoffset: z-offset

    :param upnull: toggle upper x-point

    :param lonull: toggle lower x-point

    :param npts: int
        number of points (per quadrant)

    :param doPlot: plot boundary shape construction

    :param newsq: A 4 element array, into which the new squareness values are stored

    :return: tuple with arrays of r,z,zref

    >> boundaryShape(a=0.608,eps=0.374,kapu=1.920,kapl=1.719,delu=0.769,dell=0.463,zetaou=-0.155,zetaiu=-0.255,zetail=-0.174,zetaol=-0.227,zoffset=0.000,upnull=False,lonull=False,doPlot=True)
    '''
    ukap = kapu
    lkap = kapl
    utri = delu
    ltri = dell
    uosq = zetaou
    uisq = zetaiu
    lisq = zetail
    losq = zetaol
    amin = a
    newsq[0:4] = [uosq, uisq, lisq, losq]
    #if is_int(npts):
    #    ang = np.linspace(0, 2 * np.pi, (npts * 4 + 1))
    #else:
    #    ang = npts
    ang = np.linspace(0, 2 * np.pi, (npts * 4 + 1))

    i1 = np.where((ang >= 0 * np.pi / 2.0) & (ang < 1 * np.pi / 2.0))
    i2 = np.where((ang >= 1 * np.pi / 2.0) & (ang < 2 * np.pi / 2.0))
    i3 = np.where((ang >= 2 * np.pi / 2.0) & (ang < 3 * np.pi / 2.0))
    i4 = np.where((ang >= 3 * np.pi / 2.0) & (ang < 4 * np.pi / 2.0))

    ang1 = ang[i1]
    ang2 = ang[i2]
    ang3 = ang[i3]
    ang4 = ang[i4]

    rsr2 = 1.0 / np.sqrt(2.0)
    cc = 1.0 - rsr2

    if uosq < -1.0 * rsr2:
        uosq = -1.0 * rsr2

    if uisq < -1.0 * rsr2:
        uisq = -1.0 * rsr2

    if lisq < -1.0 * rsr2:
        lisq = -1.0 * rsr2

    if losq < -1.0 * rsr2:
        losq = -1.0 * rsr2

    # n1=-alog(2.)/alog(uosq*cc+rsr2)
    n1 = -np.log(2.0) / np.log(uosq * cc + rsr2)
    # r1=amin*(1./eps-utri)+amin*(1.+utri)*cos(ang1)
    r1 = amin * (1.0 / eps - utri) + amin * (1.0 + utri) * np.cos(ang1) ** (2.0 / n1)
    # z1=zoffset+amin*ukap*sin(ang1) ^(2./n1)
    z1 = zoffset + amin * ukap * np.sin(ang1) ** (2.0 / n1)
    # z1ref=zoffset+amin*ukap*sin(ang1)
    z1ref = zoffset + amin * ukap * np.sin(ang1)
    # n2=-alog(2.)/alog(uisq*cc+rsr2)
    n2 = -np.log(2.0) / np.log(uisq * cc + rsr2)
    # r2=amin*(1./eps-utri)-amin*(1.-utri)*abs(cos(ang2))
    r2 = amin * (1.0 / eps - utri) - amin * (1.0 - utri) * abs(np.cos(ang2)) ** (2.0 / n2)
    # z2=zoffset+amin*ukap*sin(ang2) ^(2./n2)
    z2 = zoffset + amin * ukap * np.sin(ang2) ** (2.0 / n2)
    # z2ref=zoffset+amin*ukap*sin(ang2)
    z2ref = zoffset + amin * ukap * np.sin(ang2)

    if upnull:

        ##        f=findgen(99)/100.+0.01
        f = np.linspace(0.01, 1, 100)
        ##        n=findgen(99)/25.+1.04
        n = np.linspace(1.04, 5, 100)
        ##
        ##        h1=1.-(1.-uosq)*cc
        h1 = 1.0 - (1.0 - uosq) * cc
        ##        h2=1.-(1.-uisq)*cc
        h2 = 1.0 - (1.0 - uisq) * cc
        ##
        ##        a1=amin*(1.+utri)
        a1 = amin * (1.0 + utri)
        ##        a2=amin*(1.-utri)
        a2 = amin * (1.0 - utri)
        ##        b=amin*ukap
        b = amin * ukap
        ##
        ##        if (utri ge 0.) then c1=utri-1. else c1=-1./(1.+utri)
        c1 = utri - 1 if utri >= 0 else -1.0 / (1.0 + utri)
        ##        if (utri ge 0.) then c2=-1./(1.-utri) else c2=-1.*(1.+utri)
        c2 = -1.0 / (1.0 - utri) if utri >= 0.0 else -1.0 * (1.0 + utri)
        ##
        ##        y1q1=fltarr(99,99)
        y1q1 = np.zeros((100, 100))
        ##        y2q1=fltarr(99,99)
        y2q1 = np.zeros((100, 100))
        ##        y1q2=fltarr(99,99)
        y1q2 = np.zeros((100, 100))
        ##        y2q2=fltarr(99,99)
        y2q2 = np.zeros((100, 100))
        ##
        ##        for i=0,98 do begin
        for i in range(y1q1.shape[0]):
            ##        for j=0,98 do begin
            for j in range(y1q1.shape[1]):
                ##                y1q1(i,j)=(f(i)+h1*(1.-f(i))) ^n(j)+(1.-f(i) ^n(j))*h1 ^n(j)-1.
                y1q1[j, i] = (f[i] + h1 * (1.0 - f[i])) ** n[j] + (1.0 - f[i] ** n[j]) * h1 ** n[j] - 1.0
                ##                y2q1(i,j)=f(i) ^(n(j)-1.)*(f(i)*(c1+b/a1)-b/a1)-c1
                y2q1[j, i] = f[i] ** (n[j] - 1.0) * (f[i] * (c1 + b / a1) - b / a1) - c1
                ##                y1q2(i,j)=(f(i)+h2*(1.-f(i))) ^n(j)+(1.-f(i) ^n(j))*h2^n(j)-1.
                y1q2[j, i] = (f[i] + h2 * (1.0 - f[i])) ** n[j] + (1.0 - f[i] ** n[j]) * h2 ** n[j] - 1.0
                ##                y2q2(i,j)=f(i) ^(n(j)-1.)*(f(i)*(c2+b/a2)-b/a2)-c2
                y2q2[j, i] = f[i] ** (n[j] - 1.0) * (f[i] * (c2 + b / a2) - b / a2) - c2
        ##        endfor
        ##        endfor
        ##        contour,y1q1,f,n,/overplot,level=[0.0],path_xy=xy1q1,path_info=info,closed=0,/path_data_coords
        xy1q1 = contourPaths(f, n, y1q1, [0], remove_boundary_points=True)[0][0].vertices
        ##        contour,y2q1,f,n,/overplot,level=[0.0],path_xy=xy2q1,path_info=info,closed=0,/path_data_coords
        xy2q1 = contourPaths(f, n, y2q1, [0], remove_boundary_points=True)[0][0].vertices
        ##        contour,y1q2,f,n,/overplot,level=[0.0],path_xy=xy1q2,path_info=info,closed=0,/path_data_coords
        xy1q2 = contourPaths(f, n, y1q2, [0], remove_boundary_points=True)[0][0].vertices
        ##        contour,y2q2,f,n,/overplot,level=[0.0],path_xy=xy2q2,path_info=info,closed=0,/path_data_coords
        xy2q2 = contourPaths(f, n, y2q2, [0], remove_boundary_points=True)[0][0].vertices
        ##
        ##        y1q1sol=interpol(xy1q1(1,*),xy1q1(0,*),f)
        y1q1sol = np.interp(f, xy1q1[:, 0], xy1q1[:, 1])
        ##        y2q1sol=interpol(xy2q1(1,*),xy2q1(0,*),f)
        y2q1sol = np.interp(f, xy2q1[:, 0], xy2q1[:, 1])
        ##        y1q2sol=interpol(xy1q2(1,*),xy1q2(0,*),f)
        y1q2sol = np.interp(f, xy1q2[:, 0], xy1q2[:, 1])
        ##        y2q2sol=interpol(xy2q2(1,*),xy2q2(0,*),f)
        y2q2sol = np.interp(f, xy2q2[:, 0], xy2q2[:, 1])

        ##        maxdiffq1=max(y1q1sol-y2q1sol)
        maxdiffq1 = max(y1q1sol - y2q1sol)
        ##        mindiffq1=min(y1q1sol-y2q1sol)
        mindiffq1 = min(y1q1sol - y2q1sol)
        ##        maxdiffq2=max(y1q2sol-y2q2sol)
        maxdiffq2 = max(y1q2sol - y2q2sol)
        ##        mindiffq2=min(y1q2sol-y2q2sol)
        mindiffq2 = min(y1q2sol - y2q2sol)
        ##
        ##        if (maxdiffq1/mindiffq1 lt 0.) then begin
        if maxdiffq1 / mindiffq1 < 0.0:
            ##        y12q1sol=min(abs(y1q1sol-y2q1sol),imin)
            y12q1sol = min(abs(y1q1sol - y2q1sol))
            imin = np.argmin(abs(y1q1sol - y2q1sol))
            ##        fsolq1=f(imin)
            fsolq1 = f[imin]
            ##        nsolq1=y1q1sol(imin)
            nsolq1 = y1q1sol[imin]
            ##        gsolq1=(1.-fsolq1^nsolq1)^(1./nsolq1)
            gsolq1 = (1.0 - fsolq1**nsolq1) ** (1.0 / nsolq1)
        ##        endif else begin
        else:
            ##        if (maxdiffq1 gt 0.) then begin
            if maxdiffq1 > 0:
                ##          y1new=(f[94]+h1*(1.-f[94]))^y2q1sol[94]+(1.-f[94]^y2q1sol[94])*h1^y2q1sol[94]-1.
                y1new = (f[94] + h1 * (1.0 - f[94])) ** y2q1sol[94] + (1.0 - f[94] ** y2q1sol[94]) * h1 ** y2q1sol[94] - 1.0
                ##          y2new=f[94]^(y2q1sol[94]-1.)*(f[94]*(c1+b/a1)-b/a1)-c1
                y2new = f[94] ** (y2q1sol[94] - 1.0) * (f[94] * (c1 + b / a1) - b / a1) - c1
                ##          while (y1new gt y2new) do begin
                while y1new > y2new:
                    ##            h1=h1-0.01
                    h1 = h1 - 0.01
                    ##            y1new=(f[94]+h1*(1.-f[94]))^y2q1sol[94]+(1.-f[94]^y2q1sol[94])*h1^y2q1sol[94]-1.
                    y1new = (f[94] + h1 * (1.0 - f[94])) ** y2q1sol[94] + (1.0 - f[94] ** y2q1sol[94]) * h1 ** y2q1sol[94] - 1.0
                ##          endwhile
                ##          fsolq1=f[94]
                fsolq1 = f[94]
                ##          nsolq1=y2q1sol[94]
                nsolq1 = y2q1sol[94]
                ##          gsolq1=(1.-fsolq1^nsolq1)^(1./nsolq1)
                gsolq1 = (1.0 - fsolq1**nsolq1) ** (1.0 / nsolq1)
            ##        endif else begin
            else:
                ##          y1new=(f[4]+h1*(1.-f[4]))^y2q1sol[4]+(1.-f[4]^y2q1sol[4])*h1^y2q1sol[4]-1.
                y1new = (f[4] + h1 * (1.0 - f[4])) ** y2q1sol[4] + (1.0 - f[4] ** y2q1sol[4]) * h1 ** y2q1sol[4] - 1.0
                ##          y2new=f[4]^(y2q1sol[4]-1.)*(f[4]*(c1+b/a1)-b/a1)-c1
                y2new = f[4] ** (y2q1sol[4] - 1.0) * (f[4] * (c1 + b / a1) - b / a1) - c1
                ##          while (y1new lt y2new) do begin
                while y1new < y2new:
                    ##            h1=h1+0.01
                    h1 = h1 + 0.01
                    ##            y1new=(f[4]+h1*(1.-f[4]))^y2q1sol[4]+(1.-f[4]^y2q1sol[4])*h1^y2q1sol[4]-1.
                    y1new = (f[4] + h1 * (1.0 - f[4])) ** y2q1sol[4] + (1.0 - f[4] ** y2q1sol[4]) * h1 ** y2q1sol[4] - 1.0
                ##          endwhile
                ##          fsolq1=f[4]
                fsolq1 = f[4]
                ##          nsolq1=y2q1sol[4]
                nsolq1 = y2q1sol[4]
                ##          gsolq1=(1.-fsolq1^nsolq1)^(1./nsolq1)
                gsolq1 = (1.0 - fsolq1**nsolq1) ** (1.0 / nsolq1)
            ##        endelse
            ##        sqnew1=1.-(1.-h1)/cc
            sqnew1 = 1.0 - (1.0 - h1) / cc
            newsq[0] = sqnew1
        ##        endelse
        ##
        ##        alpha1=a1/(1.-fsolq1)
        alpha1 = a1 / (1.0 - fsolq1)
        ##        beta1=b/gsolq1
        beta1 = b / gsolq1
        ##
        ##        y1=beta1*(1.-((r1-amin*(1./eps+1.))/alpha1+1.)^nsolq1)^(1./nsolq1)
        y1 = beta1 * (1.0 - ((r1 - amin * (1.0 / eps + 1.0)) / alpha1 + 1.0) ** nsolq1) ** (1.0 / nsolq1)
        ##        z1=y1+zoffset
        z1 = y1 + zoffset
        ##
        ##        if (maxdiffq2/mindiffq2 lt 0.) then begin
        if maxdiffq2 / mindiffq2 < 0.0:
            ##        y12q2sol=min(abs(y1q2sol-y2q2sol),imin)
            y12q2sol = min(abs(y1q2sol - y2q2sol))
            imin = np.argmin(abs(y1q2sol - y2q2sol))
            ##        fsolq2=f(imin)
            fsolq2 = f[imin]
            ##        nsolq2=y1q2sol(imin)
            nsolq2 = y1q2sol[imin]
            ##        gsolq2=(1.-fsolq2^nsolq2)^(1./nsolq2)
            gsolq2 = (1.0 - fsolq2**nsolq2) ** (1.0 / nsolq2)
        ##        endif else begin
        else:
            ##        if (maxdiffq2 gt 0.) then begin
            if maxdiffq2 > 0.0:
                ##          y1new=(f[94]+h2*(1.-f[94]))^y2q2sol[94]+(1.-f[94]^y2q2sol[94])*h2^y2q2sol[94]-1.
                y1new = (f[94] + h2 * (1.0 - f[94])) ** y2q2sol[94] + (1.0 - f[94] ** y2q2sol[94]) * h2 ** y2q2sol[94] - 1.0
                ##          y2new=f[94]^(y2q2sol[94]-1.)*(f[94]*(c2+b/a2)-b/a2)-c2
                y2new = f[94] ** (y2q2sol[94] - 1.0) * (f[94] * (c2 + b / a2) - b / a2) - c2
                ##          while (y1new gt y2new) do begin
                while y1new > y2new:
                    ##            h2=h2-0.01
                    h2 = h2 - 0.01
                    ##            y1new=(f[94]+h2*(1.-f[94]))^y2q2sol[94]+(1.-f[94]^y2q2sol[94])*h2^y2q2sol[94]-1.
                    y1new = (f[94] + h2 * (1.0 - f[94])) ** y2q2sol[94] + (1.0 - f[94] ** y2q2sol[94]) * h2 ** y2q2sol[94] - 1.0
                ##          endwhile
                ##          fsolq2=f[94]
                fsolq2 = f[94]
                ##          nsolq2=y2q2sol[94]
                nsolq2 = y2q2sol[94]
                ##          gsolq2=(1.-fsolq2^nsolq2)^(1./nsolq2)
                gsolq2 = (1.0 - fsolq2**nsolq2) ** (1.0 / nsolq2)
            ##          endif else begin
            else:
                ##          y1new=(f[4]+h2*(1.-f[4]))^y2q2sol[4]+(1.-f[4]^y2q2sol[4])*h2^y2q2sol[4]-1.
                y1new = (f[4] + h2 * (1.0 - f[4])) ** y2q2sol[4] + (1.0 - f[4] ** y2q2sol[4]) * h2 ** y2q2sol[4] - 1.0
                ##          y2new=f[4]^(y2q2sol[4]-1.)*(f[4]*(c2+b/a2)-b/a2)-c2
                y2new = f[4] ** (y2q2sol[4] - 1.0) * (f[4] * (c2 + b / a2) - b / a2) - c2
                ##          while (y1new lt y2new) do begin
                while y1new < y2new:
                    ##            h2=h2+0.01
                    h2 = h2 + 0.01
                    ##            y1new=(f[4]+h2*(1.-f[4]))^y2q2sol[4]+(1.-f[4]^y2q2sol[4])*h2^y2q2sol[4]-1.
                    y1new = (f[4] + h2 * (1.0 - f[4])) ** y2q2sol[4] + (1.0 - f[4] ** y2q2sol[4]) * h2 ** y2q2sol[4] - 1.0
                ##          endwhile
                ##          fsolq2=f[4]
                fsolq2 = f[4]
                ##          nsolq2=y2q2sol[4]
                nsolq2 = y2q2sol[4]
                ##          gsolq2=(1.-fsolq2^nsolq2)^(1./nsolq2)
                gsolq2 = (1.0 - fsolq2**nsolq2) ** (1.0 / nsolq2)
            ##          endelse
            ##        sqnew2=1.-(1.-h2)/cc
            sqnew2 = 1.0 - (1.0 - h2) / cc
            newsq[1] = sqnew2
        ##        endelse

        ##
        ##        alpha2=a2/(1.-fsolq2)
        alpha2 = a2 / (1.0 - fsolq2)
        ##        beta2=b/gsolq2
        beta2 = b / gsolq2
        ##
        ##        y2=beta2*(1.-(1.+(amin*(1./eps-1.)-r2)/alpha2)^nsolq2)^(1./nsolq2)
        y2 = beta2 * (1.0 - (1.0 + (amin * (1.0 / eps - 1.0) - r2) / alpha2) ** nsolq2) ** (1.0 / nsolq2)
        ##        z2=y2+zoffset
        z2 = y2 + zoffset
    ##    endif

    # n3=-alog(2.)/alog(lisq*cc+rsr2)
    n3 = -np.log(2.0) / np.log(lisq * cc + rsr2)
    # r3=amin*(1./eps-ltri)-amin*(1.-ltri)*abs(cos(ang(180:269)))
    r3 = amin * (1.0 / eps - ltri) - amin * (1.0 - ltri) * abs(np.cos(ang3)) ** (2.0 / n3)
    # z3=zoffset-amin*lkap*abs(sin(ang(180:269))) ^(2./n3)
    z3 = zoffset - amin * lkap * abs(np.sin(ang3)) ** (2.0 / n3)
    # z3ref=zoffset+amin*lkap*sin(ang(180:269))
    z3ref = zoffset + amin * lkap * np.sin(ang3)
    # n4=-alog(2.)/alog(losq*cc+rsr2)
    n4 = -np.log(2.0) / np.log(losq * cc + rsr2)
    # r4=amin*(1./eps-ltri)+amin*(1.+ltri)*cos(ang(270:359))
    r4 = amin * (1.0 / eps - ltri) + amin * (1.0 + ltri) * abs(np.cos(ang4)) ** (2.0 / n4)
    # z4=zoffset-amin*lkap*abs(sin(ang(270:359))) ^(2./n4)
    z4 = zoffset - amin * lkap * abs(np.sin(ang4)) ** (2.0 / n4)
    # z4ref=zoffset+amin*lkap*sin(ang(270:359))
    z4ref = zoffset + amin * lkap * np.sin(ang4)

    if lonull:
        # f=findgen(99)/100.+0.01
        f = np.arange(99) / 100.0 + 0.01
        # n=findgen(99)/25.+1.04
        n = np.arange(99) / 25.0 + 1.04

        # h4=1.-(1.-losq)*cc
        h4 = 1.0 - (1.0 - losq) * cc
        # h3=1.-(1.-lisq)*cc
        h3 = 1.0 - (1.0 - lisq) * cc

        # a4=amin*(1.+ltri)
        a4 = amin * (1.0 + ltri)
        # a3=amin*(1.-ltri)
        a3 = amin * (1.0 - ltri)
        # b=amin*lkap
        b = amin * lkap

        # if (ltri ge 0.) then c4=ltri-1. else c4=-1./(1.+ltri)
        c4 = ltri - 1.0 if (ltri >= 0.0) else -1.0 / (1.0 + ltri)
        # if (ltri ge 0.) then c3=-1./(1.-ltri) else c3=-1.*(1.+ltri)
        c3 = -1.0 / (1.0 - ltri) if (ltri >= 0.0) else -1.0 * (1.0 + ltri)

        # y1q4=fltarr(99,99)
        y1q4 = np.zeros((99, 99))
        # y2q4=fltarr(99,99)
        y2q4 = np.zeros((99, 99))
        # y1q3=fltarr(99,99)
        y1q3 = np.zeros((99, 99))
        # y2q3=fltarr(99,99)
        y2q3 = np.zeros((99, 99))

        # for i=0,98 do begin
        for i in range(len(f)):
            # for j=0,98 do begin
            for j in range(len(n)):
                # y1q4(i,j)=(f(i)+h4*(1.-f(i))) ^n(j)+(1.-f(i) ^n(j))*h4 ^n(j)-1.
                y1q4[j, i] = (f[i] + h4 * (1.0 - f[i])) ** n[j] + (1.0 - f[i] ** n[j]) * h4 ** n[j] - 1.0
                # y2q4(i,j)=f(i) ^(n(j)-1.)*(f(i)*(c4+b/a4)-b/a4)-c4
                y2q4[j, i] = f[i] ** (n[j] - 1.0) * (f[i] * (c4 + b / a4) - b / a4) - c4
                # y1q3(i,j)=(f(i)+h3*(1.-f(i))) ^n(j)+(1.-f(i) ^n(j))*h3 ^n(j)-1.
                y1q3[j, i] = (f[i] + h3 * (1.0 - f[i])) ** n[j] + (1.0 - f[i] ** n[j]) * h3 ** n[j] - 1.0
                # y2q3(i,j)=f(i) ^(n(j)-1.)*(f(i)*(c3+b/a3)-b/a3)-c3
                y2q3[j, i] = f[i] ** (n[j] - 1.0) * (f[i] * (c3 + b / a3) - b / a3) - c3
        # endfor
        # endfor

        # contour,y1q4,f,n,/overplot,level=[0.0],path_xy=xy1q4,path_info=info,closed=0,/path_data_coords
        xy1q4 = contourPaths(f, n, y1q4, [0], remove_boundary_points=True)[0][0].vertices
        # contour,y2q4,f,n,/overplot,level=[0.0],path_xy=xy2q4,path_info=info,closed=0,/path_data_coords
        xy2q4 = contourPaths(f, n, y2q4, [0], remove_boundary_points=True)[0][0].vertices
        # contour,y1q3,f,n,/overplot,level=[0.0],path_xy=xy1q3,path_info=info,closed=0,/path_data_coords
        xy1q3 = contourPaths(f, n, y1q3, [0], remove_boundary_points=True)[0][0].vertices
        # contour,y2q3,f,n,/overplot,level=[0.0],path_xy=xy2q3,path_info=info,closed=0,/path_data_coords
        xy2q3 = contourPaths(f, n, y2q3, [0], remove_boundary_points=True)[0][0].vertices

        # y1q4sol=interpol(xy1q4(1,*),xy1q4(0,*),f)
        y1q4sol = np.interp(f, xy1q4[:, 0], xy1q4[:, 1])
        # y2q4sol=interpol(xy2q4(1,*),xy2q4(0,*),f)
        y2q4sol = np.interp(f, xy2q4[:, 0], xy2q4[:, 1])
        # y1q3sol=interpol(xy1q3(1,*),xy1q3(0,*),f)
        y1q3sol = np.interp(f, xy1q3[:, 0], xy1q3[:, 1])
        # y2q3sol=interpol(xy2q3(1,*),xy2q3(0,*),f)
        y2q3sol = np.interp(f, xy2q3[:, 0], xy2q3[:, 1])

        # maxdiffq4=max(y1q4sol-y2q4sol)
        maxdiffq4 = max(y1q4sol - y2q4sol)
        # mindiffq4=min(y1q4sol-y2q4sol)
        mindiffq4 = min(y1q4sol - y2q4sol)
        # maxdiffq3=max(y1q3sol-y2q3sol)
        maxdiffq3 = max(y1q3sol - y2q3sol)
        # mindiffq3=min(y1q3sol-y2q3sol)
        mindiffq3 = min(y1q3sol - y2q3sol)

        # if (maxdiffq4/mindiffq4 lt 0.) then begin
        if maxdiffq4 / mindiffq4 < 0.0:
            # y12q4sol=min(abs(y1q4sol-y2q4sol),imin)
            y12q4sol = min(abs(y1q4sol - y2q4sol))
            imin = np.argmin(abs(y1q4sol - y2q4sol))
            # fsolq4=f(imin)
            fsolq4 = f[imin]
            # nsolq4=y1q4sol(imin)
            nsolq4 = y1q4sol[imin]
            # gsolq4=(1.-fsolq4 ^nsolq4) ^(1./nsolq4)
            gsolq4 = (1.0 - fsolq4**nsolq4) ** (1.0 / nsolq4)
        # endif else begin
        else:
            # if (maxdiffq4 gt 0.) then begin
            if maxdiffq4 > 0.0:
                # y1new=(f(94)+h4*(1.-f(94))) ^y2q4sol(94)+(1.-f(94) ^y2q4sol(94))*h4 ^y2q4sol(94)-1.
                y1new = (f[94] + h4 * (1.0 - f[94])) ** y2q4sol[94] + (1.0 - f[94] ** y2q4sol[94]) * h4 ** y2q4sol[94] - 1.0
                # y2new=f(94) ^(y2q4sol(94)-1.)*(f(94)*(c4+b/a4)-b/a4)-c4
                y2new = f[94] ** (y2q4sol[94] - 1.0) * (f[94] * (c4 + b / a4) - b / a4) - c4
                # while (y1new gt y2new) do begin
                while y1new > y2new:
                    # h4=h4-0.01
                    h4 = h4 - 0.01
                    # y1new=(f(94)+h4*(1.-f(94))) ^y2q4sol(94)+(1.-f(94) ^y2q4sol(94))*h4 ^y2q4sol(94)-1.
                    y1new = (f[94] + h4 * (1.0 - f[94])) ** y2q4sol[94] + (1.0 - f[94] ** y2q4sol[94]) * h4 ** y2q4sol[94] - 1.0
                # endwhile

                # fsolq4=f(94)
                fsolq4 = f[94]
                # nsolq4=y2q4sol(94)
                nsolq4 = y2q4sol[94]
                # gsolq4=(1.-fsolq4 ^nsolq4) ^(1./nsolq4)
                gsolq4 = (1.0 - fsolq4**nsolq4) ** (1.0 / nsolq4)
            # endif else begin
            else:
                # y1new=(f(4)+h4*(1.-f(4))) ^y2q4sol(4)+(1.-f(4) ^y2q4sol(4))*h4 ^y2q4sol(4)-1.
                y1new = (f[4] + h4 * (1.0 - f[4])) ** y2q4sol[4] + (1.0 - f[4] ** y2q4sol[4]) * h4 ** y2q4sol[4] - 1.0
                # y2new=f(4) ^(y2q4sol(4)-1.)*(f(4)*(c4+b/a4)-b/a4)-c4
                y2new = f[4] ** (y2q4sol[4] - 1.0) * (f[4] * (c4 + b / a4) - b / a4) - c4
                # while (y1new lt y2new) do begin
                while y1new < y2new:
                    # h4=h4+0.01
                    h4 = h4 + 0.01
                    # y1new=(f(4)+h4*(1.-f(4))) ^y2q4sol(4)+(1.-f(4) ^y2q4sol(4))*h4 ^y2q4sol(4)-1.
                    y1new = (f[4] + h4 * (1.0 - f[4])) ** y2q4sol[4] + (1.0 - f[4] ** y2q4sol[4]) * h4 ** y2q4sol[4] - 1.0
                # endwhile

                # fsolq4=f(4)
                fsolq4 = f[4]
                # nsolq4=y2q4sol(4)
                nsolq4 = y2q4sol[4]
                # gsolq4=(1.-fsolq4 ^nsolq4) ^(1./nsolq4)
                gsolq4 = (1.0 - fsolq4**nsolq4) ** (1.0 / nsolq4)
            # endelse

            # sqnew4=1.-(1.-h4)/cc
            sqnew4 = 1.0 - (1.0 - h4) / cc
            newsq[3] = sqnew4
        # endelse

        # alpha4=a4/(1.-fsolq4)
        alpha4 = a4 / (1.0 - fsolq4)
        # beta4=b/gsolq4
        beta4 = b / gsolq4

        # y4=-1.*(beta4*(1.-((r4-amin*(1./eps+1.))/alpha4+1.) ^nsolq4) ^(1./nsolq4))
        y4 = -1.0 * (beta4 * (1.0 - ((r4 - amin * (1.0 / eps + 1.0)) / alpha4 + 1.0) ** nsolq4) ** (1.0 / nsolq4))
        # z4=y4+zoffset
        z4 = y4 + zoffset

        # if (maxdiffq3/mindiffq3 lt 0.) then begin
        if maxdiffq3 / mindiffq3 < 0.0:
            # y12q3sol=min(abs(y1q3sol-y2q3sol),imin)
            y12q3sol = min(abs(y1q3sol - y2q3sol))
            imin = np.argmin(abs(y1q3sol - y2q3sol))
            # fsolq3=f(imin)
            fsolq3 = f[imin]
            # nsolq3=y1q3sol(imin)
            nsolq3 = y1q3sol[imin]
            # gsolq3=(1.-fsolq3 ^nsolq3) ^(1./nsolq3)
            gsolq3 = (1.0 - fsolq3**nsolq3) ** (1.0 / nsolq3)
        # endif else begin
        else:
            # if (maxdiffq3 gt 0.) then begin
            if maxdiffq3 > 0.0:
                # y1new=(f(94)+h3*(1.-f(94))) ^y2q3sol(94)+(1.-f(94) ^y2q3sol(94))*h3 ^y2q3sol(94)-1.
                y1new = (f[94] + h3 * (1.0 - f[94])) ** y2q3sol[94] + (1.0 - f[94] ** y2q3sol[94]) * h3 ** y2q3sol[94] - 1.0
                # y2new=f(94) ^(y2q3sol(94)-1.)*(f(94)*(c3+b/a3)-b/a3)-c3
                y2new = f[94] ** (y2q3sol[94] - 1.0) * (f[94] * (c3 + b / a3) - b / a3) - c3
                # while (y1new gt y2new) do begin
                while y1new > y2new:
                    # h3=h3-0.01
                    h3 = h3 - 0.01
                    # y1new=(f(94)+h3*(1.-f(94))) ^y2q3sol(94)+(1.-f(94) ^y2q3sol(94))*h3 ^y2q3sol(94)-1.
                    y1new = (f[94] + h3 * (1.0 - f[94])) ** y2q3sol[94] + (1.0 - f[94] ** y2q3sol[94]) * h3 ** y2q3sol[94] - 1.0
                # endwhile

                # fsolq3=f(94)
                fsolq3 = f[94]
                # nsolq3=y2q3sol(94)
                nsolq3 = y2q3sol[94]
                # gsolq3=(1.-fsolq3 ^nsolq3) ^(1./nsolq3)
                gsolq3 = (1.0 - fsolq3**nsolq3) ** (1.0 / nsolq3)
            # endif else begin
            else:
                # y1new=(f(4)+h3*(1.-f(4))) ^y2q3sol(4)+(1.-f(4) ^y2q3sol(4))*h3 ^y2q3sol(4)-1.
                y1new = (f[4] + h3 * (1.0 - f[4])) ** y2q3sol[4] + (1.0 - f[4] ** y2q3sol[4]) * h3 ** y2q3sol[4] - 1.0
                # y2new=f(4) ^(y2q3sol(4)-1.)*(f(4)*(c3+b/a3)-b/a3)-c3
                y2new = f[4] ** (y2q3sol[4] - 1.0) * (f[4] * (c3 + b / a3) - b / a3) - c3
                # while (y1new lt y2new) do begin
                while y1new < y2new:
                    # h3=h3+0.01
                    h3 = h3 + 0.01
                    # y1new=(f(4)+h3*(1.-f(4))) ^y2q3sol(4)+(1.-f(4) ^y2q3sol(4))*h3 ^y2q3sol(4)-1.
                    y1new = (f[4] + h3 * (1.0 - f[4])) ** y2q3sol[4] + (1.0 - f[4] ** y2q3sol[4]) * h3 ** y2q3sol[4] - 1.0
                # endwhile

                # fsolq3=f(4)
                fsolq3 = f[4]
                # nsolq3=y2q3sol(4)
                nsolq3 = y2q3sol[4]
                # gsolq3=(1.-fsolq3 ^nsolq3) ^(1./nsolq3)
                gsolq3 = (1.0 - fsolq3**nsolq3) ** (1.0 / nsolq3)
            # endelse

            # sqnew3=1.-(1.-h3)/cc
            sqnew3 = 1.0 - (1.0 - h3) / cc
            newsq[2] = sqnew3
        # endelse

        # alpha3=a3/(1.-fsolq3)
        alpha3 = a3 / (1.0 - fsolq3)
        # beta3=b/gsolq3
        beta3 = b / gsolq3

        # y3=-1.*(beta3*(1.-(1.+(amin*(1./eps-1.)-r3)/alpha3) ^nsolq3) ^(1./nsolq3))
        y3 = -1.0 * (beta3 * (1.0 - (1.0 + (amin * (1.0 / eps - 1.0) - r3) / alpha3) ** nsolq3) ** (1.0 / nsolq3))
        # z3=y3+zoffset
        z3 = y3 + zoffset

    # r=             [r1,r2,r3,r4]
    r = np.hstack((r1, r2, r3, r4))
    # z=             [z1,z2,z3,z4]
    z = np.hstack((z1, z2, z3, z4))
    # zref=             [z1ref,z2ref,z3ref,z4ref]
    zref = np.hstack((z1ref, z2ref, z3ref, z4ref))

    return r, z, zref
##################################################################




def contourPaths(x, y, Z, levels, remove_boundary_points=False, smooth_factor=1):
    """
    :param x: 1D x coordinate

    :param y: 1D y coordinate

    :param Z: 2D data

    :param levels: levels to trace

    :param remove_boundary_points: remove traces at the boundary

    :param smooth_factor: smooth contours by cranking up grid resolution

    :return: list of segments
    """
    from matplotlib import path
    import scipy

    sf = int(round(smooth_factor))
    if sf > 1:
        x = scipy.ndimage.zoom(x, sf)
        y = scipy.ndimage.zoom(y, sf)
        Z = scipy.ndimage.zoom(Z, sf)

    [X, Y] = np.meshgrid(x, y)
    contour_generator = get_contour_generator(X, Y, Z, None, True, 0)

    mx = min(x)
    Mx = max(x)
    my = min(y)
    My = max(y)

    allsegs = []
    for level in levels:
        verts = contour_generator.create_contour(level)

        if isinstance(verts, tuple):
            # Matplotlib>3.4.3 and ContourPy return vertices and codes as a tuple.
            # Prior it was just vertices.
            segs = verts[0]
        else:
            segs = verts

        for i, seg in enumerate(segs):
            segs[i] = remove_adjacent_duplicates(seg)

        if not remove_boundary_points:
            segs_ = segs
        else:
            segs_ = []
            for segarray in segs:
                segarray = np.array(segarray)
                x_ = segarray[:, 0]
                y_ = segarray[:, 1]
                valid = []
                for i in range(len(x_) - 1):
                    if np.isclose(x_[i], x_[i + 1]) and (np.isclose(x_[i], Mx) or np.isclose(x_[i], mx)):
                        continue
                    if np.isclose(y_[i], y_[i + 1]) and (np.isclose(y_[i], My) or np.isclose(y_[i], my)):
                        continue
                    valid.append((x_[i], y_[i]))
                    if i == len(x_):
                        valid.append(x_[i + 1], y_[i + 1])
                if len(valid):
                    segs_.append(np.array(valid))

        segs = list(map(path.Path, segs_))
        allsegs.append(segs)
    return allsegs


# This function has been modified from the version available in OMFIT
def get_contour_generator(X, Y, Z, mask, corner_mask, nchunk):

    return contourpy.contour_generator(
        X,
        Y,
        Z,
        name='mpl2014',
        corner_mask=corner_mask,
        line_type=contourpy.LineType.SeparateCode,
        fill_type=contourpy.FillType.OuterCode,
        chunk_size=nchunk,
    )


def remove_adjacent_duplicates(x):
    """
    Remove adjacent duplicate rows in a 2D array
    :param x: original array
    :return: array with adjacent duplicates removed
    """
    y = []
    for i in range(len(x) - 1):
        if not np.allclose(x[i], x[i + 1]):
            y.append(x[i])
    y.append(x[-1])
    return np.array(y)