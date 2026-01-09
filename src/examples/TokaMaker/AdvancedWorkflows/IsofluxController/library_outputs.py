#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 10:03:47 2025

@author: mparsons
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import numpy as np
import pickle
import os




class STATES:
    
    # Initialize all parameters from input dictionary
    def __init__(self, mygs, N_its):
        
        # Declare all variables
        self.odir = None
        self.i_max = 0
        self.psin_i = 0
        self.psin_i_prev = 0
        self.psi = 0
        self.psi_n = 0
        self.R_p = np.zeros((N_its+1,))
        self.Z_p = np.zeros((N_its+1,))
        self.R_xpu = np.zeros((N_its+1,))
        self.Z_xpu = np.zeros((N_its+1,))
        self.R_xpl = np.zeros((N_its+1,))
        self.Z_xpl = np.zeros((N_its+1,))
        self.R_lim = np.zeros((N_its+1,))
        self.Z_lim = np.zeros((N_its+1,))
        self.I_coils = np.zeros([N_its+1, 13])
        
        # Modify initial values
        psi0 = mygs.get_psi(normalized = False)
        psin0 = mygs.get_psi(normalized = True)
        self.psin_i = psin0
        self.psin_i_prev = psin0
        self.psi = [psi0]
        self.psi_n = [psin0]
        
        rp, zp = mygs.get_globals()[1]
        self.R_p[0] = rp
        self.Z_p[0] = zp
        
        RZ_xp = mygs.get_xpoints()[0]
        self.R_xpu[0] = RZ_xp[0,0]
        self.Z_xpu[0] = RZ_xp[0,1]
        self.R_xpl[0] = RZ_xp[1,0]
        self.Z_xpl[0] = RZ_xp[1,1]
        
        RZ_ref = mygs.lim_point
        self.R_lim[0] = RZ_ref[0]
        self.Z_lim[0] = RZ_ref[1]
        
        coil_currents = mygs.get_coil_currents()[0]
        self.I_coils[0,:] = list(map(coil_currents.get, coil_currents.keys()))
        
        
        

    def update(self, mygs, i, nplot):
        
        self.i_max = i
        
        self.psin_i_prev = self.psin_i
        self.psin_i = mygs.get_psi(True)
        if(i % nplot == 0):
            self.psi.append(mygs.get_psi(False))
            self.psi_n.append(mygs.get_psi(True))
        
        rp, zp = mygs.get_globals()[1]
        self.R_p[i] = rp
        self.Z_p[i] = zp
        
        RZ_xp = mygs.get_xpoints()[0]
        self.R_xpu[i] = RZ_xp[0,0]
        self.Z_xpu[i] = RZ_xp[0,1]
        if(RZ_xp.shape[0] > 1):
            self.R_xpl[i] = RZ_xp[1,0]
            self.Z_xpl[i] = RZ_xp[1,1]
        else:
            self.R_xpl[i] = np.inf
            self.Z_xpl[i] = np.inf
        
        RZ_ref = mygs.lim_point
        self.R_lim[i] = RZ_ref[0]
        self.Z_lim[i] = RZ_ref[1]
        
        coil_currents = mygs.get_coil_currents()[0]
        self.I_coils[i,:] = list(map(coil_currents.get, coil_currents.keys()))
        






def save_simulation_results(inputs, mygs, PCS, STATES):
    
    
    

    # Store data in a file
    #####################################################
    
    odir = inputs['odir']
    os.makedirs(odir, exist_ok=True)
    
    fINPUTS = odir + 'INPUTS.pkl'
    with open(fINPUTS, 'wb') as file:
        pickle.dump(inputs, file)
    
    fPCS = odir + 'PCS.pkl'
    with open(fPCS, 'wb') as file:
        pickle.dump(PCS, file)
        
    fSTATES = odir + 'STATES.pkl'
    with open(fSTATES, 'wb') as file:
        pickle.dump(STATES, file)

    #####################################################
        
        
        

    # Basic plot parameters
    #####################################################
    
    i_max = STATES.i_max + 1
    
    t = PCS.t[:i_max]*1.0e3
    
    seg_labels = ['SEG01', 'SEG02', 'SEG03', 'SEG04', 'GRID1R', 'GRID1Z', 'GRID2R', 'GRID2Z']

    #####################################################
    
    
    

    # Plot vertical motion
    #####################################################

    fig, ax = plt.subplots()
    fig.set_figheight(5)

    ax.plot(t, 1.0e2 * STATES.Z_p[:i_max],'-')
    ax.set_ylabel(r'$Z_p$ / [cm]')
    ax.grid(True, which='both')

    ax.set_xlabel(r'Time / [ms]')
    plt.savefig(odir+'zpos.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################




    # Plot horizontal motion
    #####################################################

    fig, ax = plt.subplots()
    fig.set_figheight(5)

    ax.plot(t, 1.0e2 * STATES.R_p[:i_max],'-')
    ax.set_ylabel(r'$R_p$ / [cm]')
    ax.grid(True, which='both')

    ax.set_xlabel(r'Time / [ms]')
    plt.savefig(odir+'rpos.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################

    




    # Plot ISOFLUX proportional errors
    #####################################################

    fig, ax = plt.subplots()
    
    for i in range(8):
        EP = PCS.ISO_err_P[:i_max,i]
        ax.plot(t, EP/np.max(np.abs(EP)), '-')
        
    
    ax.set_ylabel('ISOFLUX Proportional Errors (Normalized)')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(seg_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'isoerrs_p.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################

    




    # Plot ISOFLUX integral errors
    #####################################################

    fig, ax = plt.subplots()
    
    for i in range(8):
        EI = PCS.ISO_err_I[:i_max,i]
        ax.plot(t, EI/np.max(np.abs(EI)), '-')
        
    
    ax.set_ylabel('ISOFLUX Integral Errors (Normalized)')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(seg_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'isoerrs_i.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################

    




    # Plot ISOFLUX derivative errors
    #####################################################

    fig, ax = plt.subplots()
    
    for i in range(8):
        ED = PCS.ISO_err_D[:i_max,i]
        ax.plot(t, ED/np.max(np.abs(ED)), '-')
        
    
    ax.set_ylabel('ISOFLUX Derivative Errors (Normalized)')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(seg_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'isoerrs_d.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################






    # Plot ISOFLUX proportional errors * gains
    #####################################################

    fig, ax = plt.subplots()
    
    for i in range(8):
        EPG = PCS.ISO_err_P[:i_max,i] * PCS.ISO_Gain_P[:i_max,i]
        ax.plot(t, EPG, '-')
        
    
    ax.set_ylabel('ISOFLUX Proportional Errors * Gains')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(seg_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'isoerrs_pg.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################






    # Plot ISOFLUX integral errors * gains
    #####################################################
    
    fig, ax = plt.subplots()
    
    for i in range(8):
        EIG = PCS.ISO_err_I[:i_max,i] * PCS.ISO_Gain_I[:i_max,i]
        ax.plot(t, EIG, '-')
        
        
    ax.set_ylabel('ISOFLUX Integral Errors * Gains')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(seg_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'isoerrs_ig.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################






    # Plot ISOFLUX derivative errors * gains
    #####################################################
    
    fig, ax = plt.subplots()
    
    for i in range(8):
        EDG = PCS.ISO_err_D[:i_max,i] * PCS.ISO_Gain_D[:i_max,i]
        ax.plot(t, EDG, '-')
        
    ax.set_ylabel('ISOFLUX Derivative Errors * Gains')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(seg_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'isoerrs_dg.png', dpi=300, bbox_inches="tight")
    plt.show()

    #####################################################






    # Plot ISOFLUX PID errors
    #####################################################
    
    fig, ax = plt.subplots()
    
    for i in range(8):
        EPG = PCS.ISO_err_P[:i_max,i] * PCS.ISO_Gain_P[:i_max,i]
        EIG = PCS.ISO_err_I[:i_max,i] * PCS.ISO_Gain_I[:i_max,i]
        EDG = PCS.ISO_err_D[:i_max,i] * PCS.ISO_Gain_D[:i_max,i]
        ax.plot(t, EPG+EIG+EDG, '-')
        
    ax.set_ylabel('ISOFLUX PID Errors / [arb]')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(seg_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'isoerrs_pid.png', dpi=300, bbox_inches="tight")
    plt.show()


    #####################################################






    # Plot coil currents
    #####################################################
    
    coil_labels = ['OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', 
                   'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL']
    
    coil_colors = plt.cm.tab20(np.linspace(0,1,20))
    
    fig, ax = plt.subplots()

    ax.plot(t, STATES.I_coils[:i_max,0],'-', color='k')
    ax.plot(t, STATES.I_coils[:i_max,1],'-', color=coil_colors[0])
    ax.plot(t, STATES.I_coils[:i_max,2],'-', color=coil_colors[2])
    ax.plot(t, STATES.I_coils[:i_max,3],'-', color=coil_colors[4])
    ax.plot(t, STATES.I_coils[:i_max,4],'-', color=coil_colors[6])
    ax.plot(t, STATES.I_coils[:i_max,5],'-', color=coil_colors[8])
    ax.plot(t, STATES.I_coils[:i_max,6],'-', color=coil_colors[10])
    ax.plot(t, STATES.I_coils[:i_max,7],'--', color=coil_colors[11])
    ax.plot(t, STATES.I_coils[:i_max,8],'--', color=coil_colors[9])
    ax.plot(t, STATES.I_coils[:i_max,9],'--', color=coil_colors[7])
    ax.plot(t, STATES.I_coils[:i_max,10],'--', color=coil_colors[5])
    ax.plot(t, STATES.I_coils[:i_max,11],'--', color=coil_colors[3])
    ax.plot(t, STATES.I_coils[:i_max,12],'--', color=coil_colors[1])
    
    ax.set_ylabel('PF Coil Currents / [kA]')
    ax.set_xlabel('Time / [ms]')
    ax.grid(True, which='both')
    plt.legend(coil_labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(odir+'PF_currents.png', dpi=300, bbox_inches="tight")
    plt.show()


    #####################################################






    # Make movie
    #####################################################
    
    print('Making movie')
    fig, ax0 = plt.subplots()
    line, = ax0.plot([], [])
    ax0.set_xlim(0, 10)
    ax0.set_ylim(-1, 1)
    
    if(np.mod(STATES.i_max,PCS.N_iso) == 0):
        N_frames = np.shape(STATES.psi_n)[0] + 1
    else:
        N_frames = np.shape(STATES.psi_n)[0] + 2
    
    def update(frame):
        
        ax0.clear()
        mygs.plot_machine(fig,ax0)
        
        if(frame == 0):
            result = STATES.psi_n[frame]
            i_frame = 0
            mygs.plot_psi(fig,ax0,psi=result,plasma_levels=[1.0],vacuum_nlevels=0,xpoint_marker=None,opoint_marker=None)
        elif(frame < N_frames - 2):
            result = STATES.psi_n[frame]
            i_frame = frame*PCS.N_iso
            mygs.plot_psi(fig,ax0,psi=result,plasma_levels=[1.0],vacuum_nlevels=0,xpoint_marker=None,opoint_marker=None)
        elif(frame == N_frames - 2):
            result = STATES.psin_i_prev
            i_frame = i_max - 2
            mygs.plot_psi(fig,ax0,psi=result,plasma_levels=[1.0],vacuum_nlevels=0,xpoint_marker=None,opoint_marker=None)
        else:
            result = STATES.psin_i
            i_frame = i_max - 1
            mygs.plot_psi(fig,ax0,psi=result,plasma_levels=[1.0],vacuum_nlevels=0,xpoint_marker=None,opoint_marker=None)
            
        
        plt.plot(STATES.R_p[i_frame], STATES.Z_p[i_frame], '+', color='k')
        
        RZ_cp = PCS.get_control_points_i(i_frame)
        for i in range(np.shape(RZ_cp)[0]):
            if(PCS.ISO_Gain_P[i_frame,i] != 0 or PCS.ISO_Gain_I[i_frame,i] != 0 or PCS.ISO_Gain_D[i_frame,i] != 0):
                plt.plot(RZ_cp[i,0], RZ_cp[i,1], '.')
                
        RZ_xp = PCS.get_xpt_targets_i(i_frame)
        for i in range(np.shape(RZ_xp)[0]):
            if(PCS.ISO_Gain_P[i_frame,4+2*i] != 0 or PCS.ISO_Gain_I[i_frame,4+2*i] != 0 or PCS.ISO_Gain_D[i_frame,4+2*i] != 0):
                plt.plot(RZ_xp[i,0], RZ_xp[i,1], '.')
                plt.plot(RZ_xp[i,0], RZ_xp[i,1], '.')
        
        plt.plot(STATES.R_xpu[i_frame], STATES.Z_xpu[i_frame], 'o', fillstyle='none', color='k')
        plt.plot(STATES.R_xpl[i_frame], STATES.Z_xpl[i_frame], 'o', fillstyle='none', color='k')

        g1_R0 = PCS.ISO_GRID1_BOUNDS[0]
        g1_Z0 = PCS.ISO_GRID1_BOUNDS[2]
        g1_dR = PCS.ISO_GRID1_BOUNDS[1] - PCS.ISO_GRID1_BOUNDS[0]
        g1_dZ = PCS.ISO_GRID1_BOUNDS[3] - PCS.ISO_GRID1_BOUNDS[2]
        grid1 = Rectangle((g1_R0, g1_Z0), g1_dR, g1_dZ, alpha=0.1)
        ax0.add_patch(grid1)
        
        g2_R0 = PCS.ISO_GRID2_BOUNDS[0]
        g2_Z0 = PCS.ISO_GRID2_BOUNDS[2]
        g2_dR = PCS.ISO_GRID2_BOUNDS[1] - PCS.ISO_GRID2_BOUNDS[0]
        g2_dZ = PCS.ISO_GRID2_BOUNDS[3] - PCS.ISO_GRID2_BOUNDS[2]
        grid2 = Rectangle((g2_R0, g2_Z0), g2_dR, g2_dZ, alpha=0.1)
        ax0.add_patch(grid2)

        ax0.set_axis_off()
        ax0.set_title(" t / [ms] = {:.2f}".format(t[i_frame]))
        
        return line,
    
    ani = animation.FuncAnimation(fig, update, frames=N_frames, interval=100)
    ani.save(odir + 'evolution.mp4') # Save as mp4, requires ffmpeg
    print('Post Processing Complete')
    #####################################################