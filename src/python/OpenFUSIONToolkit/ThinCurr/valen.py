import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.linalg import inv 



class VALENSystem:
    '''! Continuous-time valen system from Battey equation (9)
    '''

    def __init__(self, tw_torus, tw_mode, mode_drive):
        '''! Initialize a valen equation 9 model

        @param Lw  wall-wall inductance matrix
        @param Lwc  wall-coil inductance matrix
        @param Lwd  wall-plasma inductance matrix
        @param Lcw  coil-wall inductance matrix
        @param Lc  coil-coil inductance matrix
        @param Lcd  coil-plasma inductance matrix
        @param Ldw  plasma-wall inductance matrix
        @param Ldc  plasma-coil inductance matrix
        @param Ld  plasma inductance matrix
        @param Rw Wall resistance matrix
        @param Rc coil-coil resistance matrix
        @param Rp Plasma-driver resistance matrix

        @param Lw_ef  Effective reluctance modified Lw
        @param Lwc_ef  Effective reluctance modified Lwc
        @param Lwd_ef  Effective reluctance modified Lwd
        @param Lcw_ef  Effective reluctance modified Lcw
        @param Lc_ef  Effective reluctance modified Lc
        @param Lcd_ef  Effective reluctance modified Lcd
        @param Ldw_ef  Effective reluctance modified Ldw
        @param Ldc_ef  Effective reluctance modified Ldc
        @param Ld_ef  Effective reluctance modified Ld
       
        
        '''

        tw_mode.compute_Mcoil()
        tw_mode.compute_Lmat()
        tw_torus.compute_Mcoil()
        tw_torus.compute_Lmat()
        
        n_vcoils = tw_torus.n_vcoils
        torus_start= tw_torus.Lmat.shape[0] - n_vcoils
        mode_start = tw_mode.Lmat.shape[0] - n_vcoils 
        
        self.Lc = tw_torus.Lmat[torus_start:, torus_start:]

        self.Lw = tw_torus.Lmat[:torus_start,:torus_start]
        self.Mwc = tw_torus.Lmat[:torus_start,torus_start:]
        self.Mcw = self.Mwc.T

        Ld_full = tw_mode.Lmat[:mode_start,:mode_start]
        Mdc_full = tw_mode.Lmat[:mode_start,mode_start:]

        Mdw_full = tw_mode.cross_coupling(tw_torus)
        self.Mdw = mode_drive @ Mdw_full[:mode_start,:torus_start]
        self.Mwd = self.Mdw.T

        '''
        Mdw_full = tw_mode.cross_eval(tw_torus, mode_drive)
        self.Mdw = Mdw_full[:,:torus_start]
        self.Mwd = self.Mdw.T '''
        
        # mode_drive = mode_drive[:,:mode_start]
        self.Ld = mode_drive @ Ld_full @ mode_drive.T
        self.Mdc = mode_drive @ Mdc_full
        self.Mcd = self.Mdc.T
        
    
        tw_torus.compute_Rmat(copy_out=True)
        tw_mode.compute_Rmat(copy_out=True)        

        Rc = tw_torus.Rmat[torus_start:, torus_start:]
        Rw = tw_torus.Rmat[:torus_start,:torus_start]
        Rd_full = tw_mode.Rmat[:mode_start,:mode_start]
        Rd = mode_drive @ Rd_full @ mode_drive.T

        self.R = scipy.linalg.block_diag(Rw, Rc, Rd) 

    

    def compute_L_ef(self,s,a):
        
        '''
        @param a Boozer torque parameter 
        @param s Boozer stability parameter
        @param P Permiability matrix
        @param rho Reluctance matrix 
        '''

        P = (-1/(s+a*1j))*np.identity(2)
        rho = np.linalg.inv(self.Ld)@(P-1)
        
        
        Lw_ef = self.Lw + self.Mwd @ rho @ self.Mdw
        Lwd_ef = self.Mwd + self.Mwd @ rho @ self.Ld
        Ldw_ef = self.Mdw + self.Ld @ rho @ self.Mdw
        Ld_ef = self.Ld + self.Ld @ rho @ self.Ld

        Lwc_ef = self.Mwc + self.Mwd @ rho @ self.Mdc
        Ldc_ef = self.Mdc + self.Ld @ rho @ self.Mdc
        Lcw_ef = Lwc_ef.T 
        '''using my version of Lc, assuming the Lc equation in Battey is a typo!!!'''
        Lc_ef = self.Lc + self.Mcd @ rho @ self.Mdc
        Lcd_ef = Ldc_ef.T 


        self.L_ef = np.block([
            [Lw_ef, Lwc_ef, Lwd_ef],
            [Lcw_ef, Lc_ef, Lcd_ef],
            [Ldw_ef, Ldc_ef, Ld_ef]
        ])

    def eigenvalues(self,k=20):
        '''! Compute eigenvalues of the homogeneous equation 9 system

            @param eig_vals eigenvalues of equation 9 system 
            @param eig_vecs eigenvectors of equation 9 system 
        '''

        self.eig_vals, self.eig_vecs = scipy.sparse.linalg.eigs(inv(self.L_ef) @ -self.R, k, which = 'LR')
        return self.eig_vals, self.eig_vecs
    

    def plot(self):
        plt.figure()
        plt.scatter(self.eig_vals.real, self.eig_vals.imag)
        plt.axhline(0.0, color="k", linewidth=0.8)
        plt.axvline(0.0, color="k", linewidth=0.8)
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        plt.title("Eigenvalues")
        plt.grid(True)
        plt.show()
    


    
    