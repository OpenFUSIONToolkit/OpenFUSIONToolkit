import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.linalg import inv 



class RWMmodel:
    '''! Continuous-time RWM model 
    '''

    def __init__(self, tw_torus, tw_mode, mode_drive):
        '''! Calulate inductance and resistance matrices 

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

        self.Ld = mode_drive @ Ld_full @ mode_drive.T
        self.Mdc = mode_drive @ Mdc_full
        self.Mcd = self.Mdc.T
        
    
        tw_torus.compute_Rmat(copy_out=True)
        tw_mode.compute_Rmat(copy_out=True)        

        self.Rc = tw_torus.Rmat[torus_start:, torus_start:]
        self.Rw = tw_torus.Rmat[:torus_start,:torus_start]
        Rd_full = tw_mode.Rmat[:mode_start,:mode_start]
        self.Rd = mode_drive @ Rd_full @ mode_drive.T

        self.R = scipy.linalg.block_diag(self.Rw, self.Rc, self.Rd) 

    

    def compute_L_ef(self,s,a):
        
        '''! Compute effective inductance matrices with Boozer s and a parameters 

        @param a Boozer torque parameter 
        @param s Boozer stability parameter
        @param P Permiability matrix
        @param rho Reluctance matrix 

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

        P = (-1/(s+a*1j))*np.identity(2)
        self.rho = np.linalg.inv(self.Ld)@(P-np.identity(2))
        
        
        self.Lw_ef = self.Lw + self.Mwd @ self.rho @ self.Mdw
        self.Lwd_ef = self.Mwd + self.Mwd @ self.rho @ self.Ld
        self.Ldw_ef = self.Mdw + self.Ld @ self.rho @ self.Mdw
        self.Ld_ef = self.Ld + self.Ld @ self.rho @ self.Ld

        self.Lwc_ef = self.Mwc + self.Mwd @ self.rho @ self.Mdc
        self.Ldc_ef = self.Mdc + self.Ld @ self.rho @ self.Mdc
        self.Lcw_ef = self.Lwc_ef.T 
        self.Lc_ef = self.Lc + self.Mcd @ self.rho @ self.Mdc
        self.Lcd_ef = self.Ldc_ef.T 

        self.L_ef = np.block([
            [self.Lw_ef, self.Lwc_ef, self.Lwd_ef],
            [self.Lcw_ef, self.Lc_ef, self.Lcd_ef],
            [self.Ldw_ef, self.Ldc_ef, self.Ld_ef]
        ])


    def eigenvalues(self,k=20):
        '''! Compute eigenvalues of dI/dt = -L^-1R system 

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
    


    
    