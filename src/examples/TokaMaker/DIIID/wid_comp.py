#!/usr/bin/env Python
#
# Copyright (C) 2012 Thomas H. Osborne (osborne@fusion.gat.com)
#
# This program is part of the pyD3D package
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"""
---------------------------------------
Module wid_comp : Compare EPED1 Calculation of Width with Measure width
---------------------------------------
Methods:
-------
  wc(...):  
     Calculate widths in flux space
-------
  pl(...): Plot widths
-------
  plw(...): Calculat and plot comparison between measured pedestal widths and EPED1 prediction
-------
"""

import sys
__version__ = "rel_4-70-9"

def wc(
        shot,             # Shot number
        da = [],          # List of Dalpha signals to also read for comparison
        smooth = 0,       # Smooth results by this many ms
        ceped = 0.089,    # EPED1 scaling factor. w=0.076*(2*betae_pol_ped)**0.5
        tmin = 0.0,       # Start time in ms
        tmax = 5999.,     # End time in ms
        efittree='JT_TS', # EFIT tree name or runtag for mapping to psin
        use_tewid = False,# If True use Te width rather than (Tewid+newid)/2
        use_pewid = False,# If True use scale_pewid*(pe width)rather than (Tewid+newid)/2
        scale_pewid = 1.3,# If use_pewid=True use scale_pewid*(pe width)rather than (Tewid+newid)/2
        x_coord = 'psin',   # Plot resuts in this coordinate:
                            # 'psin' = normalized poloidal flux
                            # 'zts' = z along the TS chord in m
                            # 'rmid' = R along the midplane in m
        betap_from = 'pe',  # Compute betapped from either 'pe' or 'ne*te'
):
    """
    wc(
    # Compute widths in flux space
        shot,             # Shot number
        da = [],          # List of Dalpha signals to also read for comparison
        smooth = 0,       # Smooth results by this many ms
        ceped = 0.089,    # EPED1 scaling factor. w=0.076*(2*betae_pol_ped)**0.5
        tmin = 0.0,       # Start time in ms
        tmax = 5999.,     # End time in ms
        efittree='JT_TS', # EFIT tree name or runtag for mapping to psin
        use_tewid = False,# If True use Te width rather than (Tewid+newid)/2
        use_pewid = False,# If True use scale_pewid*(pe width)rather than (Tewid+newid)/2
        scale_pewid = 1.3,# If use_pewid=True use scale_pewid*(pe width)rather than (Tewid+newid)/2
        x_coord = 'psin',   # Plot resuts in this coordinate:
                            # 'psin' = normalized poloidal flux
                            # 'zts' = z along the TS chord in m
                            # 'rmid' = R along the midplane in m
        betap_from = 'pe',  # Compute betapped from either 'pe' or 'ne*te'
    )
    Returns a dictionary of widths and other data
    da = [], # List of Dalpha signals to also read
    """
    mu0 = 4.e-7*numpy.pi
    e = 1.6022e-19

    d = {}

    gdat = efit_eqdsk.get_gdat(
        shot,efit=efittree,
        gnames=['psirz','ssimag','ssibry','gtime'],quiet=1)
    adat = efit_eqdsk.get_adat(
        shot,efit=efittree,
        anames=['zuperts','z0','bpolav','rmidout'],quiet=1)

    d['tewid'] = Data('tewid',shot,quiet=1).xslice((0,tmin,tmax))
    d['tesym'] = Data('tesym',shot,quiet=1).xslice((0,tmin,tmax))
    d['newid'] = Data('newid',shot,quiet=1).xslice((0,tmin,tmax))
    d['pewid'] = Data('pewid',shot,quiet=1).xslice((0,tmin,tmax))
    d['peped'] = Data('peped',shot,quiet=1).xslice((0,tmin,tmax))
    d['teped'] = Data('teped',shot,quiet=1).xslice((0,tmin,tmax))
    d['neped'] = Data('neped',shot,quiet=1).xslice((0,tmin,tmax))
    
    if betap_from == 'pe':
        d['betapped'] = (2.e3*2.*mu0)*d['peped']/adat['bpolav']**2
        print('wid_comp.wc:', end=' ')
        print('Computing betapped using peped')
    else:
        d['betapped'] = (2.*2.*mu0*e)*d['teped']*d['neped']/adat['bpolav']**2
        print('wid_comp.wc:', end=' ')
        print('Computing betapped using neped*Teped')
    print('wid_comp.wc:', end=' ')
    print('Using C_EPED1 = %6.4f'%(ceped))
    d['weped'] = ceped*d['betapped']**0.5
    if use_tewid:
        d['wnt'] = d['tewid']
    elif use_pewid:
        d['wnt'] = scale_pewid*d['pewid']
    else:
        d['wnt'] = 0.5*(d['tewid']+d['newid'])
    d['wp']  = scale_pewid*d['pewid']
    d['ceped'] = d['wnt']/d['betapped']**0.5

    d['dpsindz'] = Data()
    d['dpsindz'].x = [ gdat['gtime'].x[0].copy() ]
    d['dpsindz'].y =   gdat['gtime'].y.copy()
    d['dpsindz'].yname = 'dpsindz'

    d['dpsindr'] = Data()
    d['dpsindr'].x = [ gdat['gtime'].x[0].copy() ]
    d['dpsindr'].y =   gdat['gtime'].y.copy()
    d['dpsindr'].yname = 'dpsindr'

    for k in range(gdat['gtime'].y.size):
        sm = gdat['ssimag'].y[k]
        sb = gdat['ssibry'].y[k]
        rts = 1.94
        zmid = adat['z0'].y[k]
        zped = adat['zuperts'].y[k]-1.0e-2
        rped = adat['rmidout'].y[k]-0.5e-2
        pn = ((gdat['psirz'][k,:,:]-sm)/(sb-sm)).spline(quiet=1)
        d['dpsindz'].y[k] = pn( rts,zped,0,1)
        d['dpsindr'].y[k] = pn(rped,zmid,1,0)

    print('wid_comp.wc:', end=' ')
    print('Average, Median dpsindz = %5.3f, %5.3f'%(
        numpy.mean(   d['dpsindz'].xslice((0,tmin,tmax)).y ),
        numpy.median( d['dpsindz'].xslice((0,tmin,tmax)).y )
        ))

    print('wid_comp.wc:', end=' ')
    print('Average, Median dpsindr = %5.3f, %5.3f'%(
        numpy.mean(   d['dpsindr'].xslice((0,tmin,tmax)).y ),
        numpy.median( d['dpsindr'].xslice((0,tmin,tmax)).y )
        ))

    if x_coord == 'psin':
        d['wnt'] *= d['dpsindz']
        d['wp']  *= d['dpsindz']
        d['ceped'] *=d['dpsindz']
        d['peprim'] = d['peped']/d['wp']
        d['peprim_eped'] = d['peped']/d['weped']
    elif x_coord == 'zts':
        d['weped'] = d['weped']/d['dpsindz']
        d['peprim'] = d['peped']/d['wp']
        d['peprim_eped'] = d['peped']/d['weped']
    elif x_coord == 'rmid':
        d['wnt']*=d['dpsindz']/d['dpsindr']
        d['wp']*=d['dpsindz']/d['dpsindr']
        d['ceped']*=d['dpsindz']/d['dpsindr']
        d['weped'] = d['weped']/d['dpsindr']
        d['peprim'] = d['peped']/d['wp']
        d['peprim_eped'] = d['peped']/d['weped']
    else:
        raise ValueError('x_coord must be one of psin,zts,rmid')
        
    for k in list(d.keys()):
        d[k].shot=shot
        if smooth > 0:
            d[k] = (d[k].newx(1.)).smooth(float(smooth),'median')
    for k in da:
        d[k] = Data(k,shot,quiet=1)

    return d

def pl(
        d,            # Dictionary of widths in flux space
        tmin=1000.,   # Start time of plot
        tmax=5999.,   # End time of plot
        wmin=0.0,     # Minimum width in plot (%psin)
        wmax=0.1,     # Maximum width in plot (%psin)
        ppmin=0.,     # Minimum pprime in plot
        ppmax=60.,    # Maximum pprime in plot
        da = [],      # List of Dalpha signals to also plot
        include_wp = False, # If True include pessure width as 1.3*prmtan_pewid
        include_gp = False, # If True plot gradp for measured and eped widths 
        include_ceped=False,# If True plot the ratio w/(2*betae_pol_ped)**0.5
        number_size = 2, # Number size in graph
        label = 'auto',# Set to None to turn off labels
        color = 'red', # Color of the experimental width in the plot
        cursor = 0,    # If 1 enter cursor mode after plotting
        ps_output = False,# If True output to color postscript file (.cps) rather than screen
):
    """
    pl(
    # Plot comparison of widths in flux space with EPED1 predictions
        d,            # Dictionary of widths in flux space
        tmin=1000.,   # Start time of plot
        tmax=5999.,   # End time of plot
        wmin=0.0,     # Minimum width in plot (%psin)
        wmax=0.1,     # Maximum width in plot (%psin)
        ppmin=0.,     # Minimum pprime in plot
        ppmax=60.,    # Maximum pprime in plot
        da = [],      # List of Dalpha signals to also plot
        include_wp = False, # If True include pessure width as 1.3*prmtan_pewid
        include_gp = False, # If True plot gradp for measured and eped widths 
        include_ceped=False,# If True plot the ratio w/(2*betae_pol_ped)**0.5
        number_size = 2, # Number size in graph
        label = 'auto',  # Set to None to turn off labels
        color = 'red', # Color of the experimental width in the plot
        cursor = 0,    # If 1 enter cursor mode after plotting
        ps_output = False,# If True output to color postscript file (.cps) rather than screen
    )
    """

    s=Screen(__name__,d)
    s.aw(xmin=tmin,xmax=tmax,iden=0)
    s.ag(ymin=wmin,ymax=wmax,number_size=number_size)
    s.ad(d['wnt'],label=label,color=color)
    if include_wp:
        s.ad(d['wp'],label=label,color='blue')
    if include_ceped:
        s.ad(d['ceped'],label=label,color='blue')
    s.ad(d['weped'],color='gray',label=label)
    for k in da:
        d[k].y = d[k].y*(0.5*wmax/d[k].y.max())
        s.ad(d[k],color='cyan',label=label)

    if include_gp:
        s.ag(ymin=ppmin,ymax=ppmax,number_size=number_size)
        s.ad(d['peprim'],color='red',line=None,symbol='filled-circle',symbol_size=4)
        s.ad(d['peprim_eped'],color='blue')

#    for k in da:s.ag(ymin=0,ymax=1,number_size=number_size);s.ad(d[k])

    if ps_output:
        print('Writing plot to PostScript file wid_comp_%i_0.cps'%(d['weped'].shot))
        s.pl('cps','wid_comp_%i_'%(d['weped'].shot),hc_match=1)
    else:
        s.pl(cursor=cursor)

def plw(
    # -----------------------------------------------------------------
    # Compute widths in flux space, and 
    # plot comparison of widths in flux space with EPED1 predictions
    # -----------------------------------------------------------------
    shot,               # Shot number
    efittree ='JT_TS',  # EFIT tree name or runtag for mapping to psin
    da = [],            # List of Dalpha signals to also plot
    use_tewid = False,  # If True use Te width rather than (Tewid+newid)/2
    use_pewid = False,  # If True use scale_pewid*(pe width)rather than (Tewid+newid)/2
    scale_pewid = 1.3,  # Scale factor for pewid
    ceped = 0.089,      # EPED1 scaling factor. w=ceped*(2*betae_pol_ped)**0.5,
                        # experimental value is 0.076, theoretical is 0.089
    x_coord = 'psin',   # Plot resuts in this coordinate:
                        # 'psin' = normalized poloidal flux
                        # 'zts' = z along the TS chord in m
                        # 'rmid' = R along the midplane in m
    betap_from = 'pe',  # Compute betapped from either 'pe' or 'ne*te'
    include_gp = False, # If True plot gradp for measured and eped widths 
    include_ceped=False,# If True plot the ratio w/(2*betae_pol_ped)**0.5
    # Plot control
    tmin = 1000.,       # Start time of plot
    tmax = 5999.,       # End time of plot
    wmin = 0.0,         # Minimum width in plot (%psin)
    wmax = 0.1,         # Maximum width in plot (%psin)
    ppmin = 0.,         # Minimum pprime in plot
    ppmax = 60.,        # Maximum pprime in plot
    smooth = 0,         # Smooth results by this many ms
    number_size = 2,    # Number size in graph
    label = 'auto',     # Set to None to turn off labels
    color = 'red',      # Color of the experimental width in the plot
    cursor = 0,         # If 1 enter cursor mode after plotting
    ps_output = False,  # If True output to color postscript file (.cps) rather than screen
    ):
    """
    plw(
    # -----------------------------------------------------------------
    # Compute widths in flux space, and 
    # plot comparison of widths in flux space with EPED1 predictions
    # -----------------------------------------------------------------
    shot,               # Shot number
    efittree ='JT_TS',  # EFIT tree name or runtag for mapping to psin
    da = [],            # List of Dalpha signals to also plot
    use_tewid = False,  # If True use Te width rather than (Tewid+newid)/2
    use_pewid = False,  # If True use scale_pewid*(pe width)rather than (Tewid+newid)/2
    scale_pewid = 1.3,  # Scale factor for pewid
    ceped = 0.089,      # EPED1 scaling factor. w=ceped*(2*betae_pol_ped)**0.5,
                        # experimental value is 0.076, theoretical is 0.089
    x_coord = 'psin',   # Plot resuts in this coordinate:
                        # 'psin' = normalized poloidal flux
                        # 'zts' = z along the TS chord in m
                        # 'rmid' = R along the midplane in m
    betap_from = 'pe',  # Compute betapped from either 'pe' or 'ne*te'
    include_gp = False, # If True plot gradp for measured and eped widths 
    include_ceped=False,# If True plot the ratio w/(2*betae_pol_ped)**0.5
    # Plot control
    tmin = 1000.,       # Start time of plot
    tmax = 5999.,       # End time of plot
    wmin = 0.0,         # Minimum width in plot (%psin)
    wmax = 0.1,         # Maximum width in plot (%psin)
    ppmin = 0.,         # Minimum pprime in plot
    ppmax = 60.,        # Maximum pprime in plot
    smooth = 0,         # Smooth results by this many ms
    number_size = 2,    # Number size in graph
    label = 'auto',     # Set to None to turn off labels
    color = 'red',      # Color of the experimental width in the plot
    cursor = 0,         # If 1 enter cursor mode after plotting
    ps_output = False,  # If True output to color postscript file (.cps) rather than screen
    )
    Returns a dictionary of widths and other data
    """
    if type(da) is not list:da=[da]
    d = wc(
        shot,da,smooth,tmin=tmin,tmax=tmax,
        ceped=ceped,efittree=efittree,
        use_tewid=use_tewid,use_pewid=use_pewid,
        scale_pewid=scale_pewid,
        x_coord = x_coord,
        betap_from = betap_from,
        )
    pl(
        d,tmin,tmax,wmin,wmax,ppmin,ppmax,da,
        number_size=number_size,
        label = label,
        color = color,
        cursor = cursor,
        ps_output = ps_output,
        include_gp=include_gp,
        include_ceped=include_ceped,
        )
    return d

################################################################################################
################################ RUN AS SCRIPT #################################################
################################################################################################
def _argtype( a ):
    if a[0] in ['(','[']:
        e = a[1:-1].split(',')
        if e[0].isalpha():return a[0]+"'"+"','".join(e)+"'"+a[-1]
        return a
    try:
        float(a)
    except:
        if a.upper()=='NONE':
            return 'None'
        elif a.upper() == 'FALSE':
            return 'False'
        elif a.upper() == 'TRUE':
            return 'True'
        else:
            return "'"+a+"'"
    else:
        return a
################################################################################################
if __name__ == '__main__':
    _methods = {
        'r': 'plw',
        }
    _descripts = {
        'r': 'Compute ped widths in flux space and plot comparison to EPED1',
        }       
    _method_order=['r']
    _doc = """
    SYNTAX:
     %s -x required_arg0_value required_arg1_value ... optional_arg0=value optional_arg1=value ...
         Where x is a switch setting run type, and one of %s. 
         Arguments must be separated by space. 
         Required argumensts must be first, in order, and without argname. 
         Optional arguments may be listed in order without the argument name or in any order
         as argname=argvalue. No space is allowed around = .
         String argument values do not need to be in quotes. None argvalue is None
         When using csh or tcsh, list type arguments must be in quotes, e.g. host_types='[hp,osf]'
      Example: %s -s 98893 tstart=2000 dt=10 nt=100 snap=jt
    """ % (sys.argv[0],_method_order,sys.argv[0])
    _doc = _doc[:] + \
        '\n    RUN CONTROL SWITCHES (For detailed help on a given run mode enter %s -hx )\n'%(sys.argv[0])
    for _k in _method_order:
        _doc = _doc[:] + "      -%s : %s\n" %(_k,_descripts[_k])
    if len(sys.argv) == 1:sys.argv = [sys.argv[0],'-h']
    args = []
    for arg in sys.argv[2:]:
        if '=' in arg:
            sarg = arg.split('=')
            args.append( sarg[0] + "=" + _argtype( sarg[1] ) )
        else:
            args.append( _argtype( arg ) )
    _meth = sys.argv[1][1:]
    print('='*100)
    print('%s: Version = %s' % (sys.argv[0],__version__))
    if 'h' in _meth:
        print('='*100)
        print(_doc)
        if len(_meth) > 1:
            _kms = _meth[ 1 - _meth.index('h') ]
            print('\n%s -%s %s\n' % ('<'*50,_kms,'>'*50))
            if _kms == 's':
                setupdb_efit(None,print_doc=1)
            elif _kms in list(_methods.keys()):
                print("  %s: %s" % ( _methods[_kms], _descripts[_kms] ))
                exec( "print (%s.__doc__)" %( _methods[_kms] ) )
            else:
                print("  RUN MODE %s IS UNKNOWN !\n" %(_kms))
    elif _meth in list(_methods.keys()):
        from data import *
        from screens import *
        import efit_eqdsk
        if 'SHOT' == _descripts[_meth][-4:]:
            _shot = int(sys.argv[2])
            print("%s: %s %i" % ( _methods[_meth], _descripts[_meth], _shot ))
        else:
            print("%s: %s" % ( _methods[_meth], _descripts[_meth] ))
        com = "%s(%s)" % ( _methods[_meth], ",".join(args) )
        print(com)
        exec(com)
    elif _meth != 'v':
        print("\n  RUN MODE %s IS UNKNOWN !\n" %(_meth))
else:
    from data import *
    from screens import *
    import efit_eqdsk
################################################################################################
##################################### THE END ##################################################
################################################################################################

