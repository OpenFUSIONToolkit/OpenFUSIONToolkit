'''! Python interface for ThinCurr thin-wall eddy current functionality

@authors Thomas Wang
@date June 2025
@ingroup doxy_oft_python
'''
import numpy as np
from scipy.__config__ import show
from sympy import N
from .sensor import Mirnov, save_sensors
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from ..io import histfile

class torus_fourier_sensor():
    '''! Interface that sets up normal magnetic field probes at a desired surface, Fourier analyze the signals, and output them for GPEC'''
    
    def __init__(self,radial_positions,axial_positions,major_radius,helicity):
        '''! Initialize the interface with surface

        @param radial_positions The R coordinates that defines the surface [ntheta]
        @param axial_positions The Z coordinates that defines the surface [ntheta]
        @param major_radius The major radius of the device
        @param helicity Sign of the plasma helicity
        '''

        if helicity not in [1,-1]:
            raise ValueError('Helicity should be either 1 or -1')
        self.major_radius = major_radius
        self.helicity = helicity
        if radial_positions[0] == radial_positions[-1] and axial_positions[0] == axial_positions[-1]:
            self.radial_positions =  radial_positions[0:-1]
            self.axial_positions = axial_positions[0:-1]
        else:
            self.radial_positions =  radial_positions
            self.axial_positions = axial_positions
        theta_values,_ = self.convert_to_polar()
        sorted_indices = np.argsort(theta_values)
        self.radial_positions = self.radial_positions[sorted_indices]
        self.axial_positions = self.axial_positions[sorted_indices]
        self.theta_list = theta_values[sorted_indices]
        self.hist_file = None

    def convert_to_polar(self):
        '''!Converts (R,Z) to (theta,r) with respect to a major radius'''
        r = np.sqrt((self.radial_positions-self.major_radius)**2+self.axial_positions**2)
        theta = np.arctan2(self.axial_positions, self.radial_positions - self.major_radius)
        theta = np.mod(theta,2*np.pi)
        return theta, r
    
    def place_normal_sensors(self,nphi=15,filename='floops.loc',ax=None):
        '''! Place normal sensors on the ThinCurr object

        @param nphi The number of poloidal planes to be probed by sensors in the toroidal direction
        @param filename The name of the file storing ThinCurr sensors
        @param ax Matplotlib axis for plotting
        '''
        
        def calculate_outward_unit_normals(R,Z,R_major):
            '''! Helper function that calculates outward-pointing unit normal vectors for a poloidal cross-section

            @param R Radial positions of the surface [ntheta]
            @param Z Axial positions of the surface [ntheta]
            @param R_major The major Radius of the surface
            @result `normals` The outward unit normal vectors on the surface [ntheta,2]
            '''

            n = len(R)
            dR = R[(np.arange(n)+1)%n]-R
            dZ = Z[(np.arange(n)+1)%n]-Z
            T = np.stack([dR,dZ],axis=1)
            
            T_magn = np.linalg.norm(T,axis=1,keepdims=True)
            T_magn[T_magn == 0] = 1e-12
            
            T_unit = T/T_magn
            
            normals = np.column_stack((T_unit[:,1],-T_unit[:,0]))

            center_vectors = np.stack([R-R_major,Z],axis=1)
            dot_products = np.einsum('ij,ij->i',normals,center_vectors)
            normals[dot_products<0] *= -1
            
            return normals

        outward_normals = calculate_outward_unit_normals(self.radial_positions,self.axial_positions,self.major_radius)

        sensors = []
        ntheta = len(self.radial_positions)
        self.ntheta = ntheta
        self.nphi = nphi
        for i, phi in enumerate(np.linspace(2*np.pi,0,nphi,endpoint=False)):
            for j in range(ntheta):
                sensors.append(Mirnov(np.r_[self.radial_positions[j]*np.cos(phi),self.radial_positions[j]*np.sin(phi),self.axial_positions[j]],np.r_[outward_normals[j,0]*np.cos(phi),outward_normals[j,0]*np.sin(phi),outward_normals[j,1]],'B_{0}_{1}'.format(i+1,j+1)))

        save_sensors(sensors,filename=filename)
        
        if ax is None:
            _,ax = plt.subplots(figsize=(10, 6))
            show_plot = True
        else:
            show_plot = False

        ax.plot(np.concatenate([self.radial_positions,self.radial_positions[:1]]),np.concatenate([self.axial_positions,self.axial_positions[:1]]))
        ax.quiver(self.radial_positions,self.axial_positions,outward_normals[:,0],outward_normals[:,1],color='red',scale=20)
        ax.set_xlabel("Radial Position (R)")
        ax.set_ylabel("Axial Position (Z)")
        ax.set_ylim((-1.25*max(self.axial_positions),1.25*max(self.axial_positions)))
        ax.set_xlim((min(self.radial_positions)-0.1,max(self.radial_positions)+0.1))
        ax.set_title("Surface with Outward Normal Vectors")
        ax.grid(True)
        ax.set_aspect('equal')
        if show_plot:
            plt.show()

    def load_histfile(self, hist_file_path='floops.hist'):
        '''! Loading histfile containing magnetic values at the surface collected by sensors

            @param hist_file_path Path to the histfile containing sensor values
        '''
        
        hist_file = histfile(hist_file_path)
        if hist_file.nfields-1 != self.ntheta*self.nphi:
            raise ValueError("The hist file provided might not be the output of the sensors defined in the current instance of the sensor_interface class.")
        else:
            self.hist_file = hist_file

    def get_B_mesh(self,t):
        '''! Extract the mesh of magnetic sensor values at a t time step

            @param t The time step at which the sensor values are extracted from sensor_file
            @result `sensor_mesh` Sensor values on a (theta,phi) mesh [ntheta,nphi]
        '''

        if self.hist_file is None:
            raise AttributeError('Probe information not available, see "load_histfile".')
        elif t >= len(self.hist_file['B_1_1']):
            raise ValueError(f'Time step is larger than the maximum time step ran in the simulation ({len(self.hist_file["B_1_1"])-1}).')
        else:
            sensor_mesh = np.zeros((self.ntheta,self.nphi))
            for i in range(self.ntheta):
                for j in range(self.nphi):
                    key = f'B_{j+1}_{i+1}'
                    sensor_mesh[i][j] = self.hist_file[key][t]
            return sensor_mesh
    
    def fft2(self,B,theta_coord='PEST',d_phi=None):
        '''! Performs a 2D Fast Fourier Transform in PEST/Hamada coordinates on the magnetic field matrix probed by sensors with proper Nyquist handling

        @param B Input matrix of shape [ntheta,nphi]
        @param theta_coord The coordinate system used for theta
        @param d_phi Hamada phase shifts [ntheta]
        @result `B_n` The Fast Fourier Transformed matrix [ntheta,nphi], `n_modes` The toroidal mode numbers [nphi], `m_modes` The poloidal mode numbers [ntheta]
        '''

        n_theta, n_phi = B.shape
        n_modes = np.round(np.fft.fftfreq(n_phi)*n_phi).astype(int)
        m_modes = np.round(np.fft.fftfreq(n_theta)*n_theta).astype(int)

        if theta_coord=='PEST':
            B_n = np.fft.fft2(B,norm='forward')
        elif theta_coord=='Hamada':
            if d_phi is not None:
                B_n_toroidal = np.fft.fft(B,axis=1,norm='forward')
                phase_matrix = np.exp(-1j*np.outer(d_phi,n_modes))
                B_n_toroidal_shifted = B_n_toroidal * phase_matrix
                B_n = np.fft.fft(B_n_toroidal_shifted,axis=0,norm='forward')
            else:
                raise ValueError('List of phase shifts in theta is required by Hamada coordinate system.')
        else:
            raise ValueError(f"Input coordinate system '{theta_coord}' is not supported. Supported systems are PEST and Hamada.")

        B_n *= 2.0
        B_n[0,0] /= 2.0
        if n_theta%2==0:
            B_n[n_theta//2,0] /= 2.0
        if n_phi%2==0:
            B_n[0,n_phi//2] /= 2.0
        if n_theta%2==0 and n_phi%2==0:
            B_n[n_theta//2,n_phi//2] /=2.0

        return B_n, n_modes, m_modes

    def ifft2(self,B_n,theta_coord='PEST',d_phi=None):
        '''! Performs an inverse 2D Fast Fourier Transform in PEST/Hamada coordinates on the transformed magnetic field matrix with proper Nyquist handling

        @param B_n Input FFT'ed matrix of shape [ntheta,nphi]
        @param theta_coord The coordinate system used for theta
        @param d_phi Hamada phase shifts [ntheta]
        @result `B_ifft` The reconstructed matrix [ntheta,nphi]
        '''

        n_theta, n_phi = B_n.shape

        B_n /= 2.0
        B_n[0,0] *= 2.0
        if n_theta%2 == 0:
            B_n[n_theta//2][0] *= 2.0
        if n_phi%2 == 0:
            B_n[0][n_phi//2] *= 2.0
        if n_theta%2 == 0 and n_phi%2 == 0:
            B_n[n_theta//2][n_phi//2] *=2.0

        if theta_coord=='PEST':
            return np.fft.ifft2(B_n,norm="forward")
        elif theta_coord=='Hamada':
            if d_phi is not None:
                B_ifft_poloidal_shifted = np.fft.ifft(B_n,axis=0,norm='forward')
                n_modes = np.round(np.fft.fftfreq(n_phi)*n_phi).astype(int)
                inverse_phase = np.exp(1j*np.outer(d_phi,n_modes))
                B_ifft_poloidal = B_ifft_poloidal_shifted*inverse_phase

                return np.fft.ifft(B_ifft_poloidal,axis=1,norm='forward')
            else:
                raise ValueError('List of phase shifts in theta is required by Hamada coordinate system.')
        else:
            raise ValueError(f"Input coordinate system '{theta_coord}' is not supported. Supported systems are PEST and Hamada.")
    
    def plot_inverse_2D_fourier_transform(self,t,fig=None,ax=None):
        '''! Use ifft2 to reconstruct the sensor signal at time step t

        @param t The time step of the sensor signals at which fft and ifft is to be performed
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        '''

        B_n_fft, _, _ = self.fft2(self.get_B_mesh(t))
        B_n_ifft = self.ifft2(B_n_fft)
        # print(f'The largest difference between the original and reconstructed sensor signals is {np.linalg.norm(abs(B_n_ifft.real-self.get_B_mesh(t).real),np.inf)}')
        phi_grid, theta_grid = np.meshgrid(np.linspace(0,2*np.pi,self.nphi,endpoint=False),self.theta_list)
        if ax is None and fig is None:
            fig,ax = plt.subplots(figsize=(10,6),constrained_layout=True)
            show_plot = True
        elif ax is not None and fig is not None:
            show_plot = False
        else:
            raise ValueError('For customized plotting, please provide both fig and ax as inputs.')
        
        ax.set_title(fr"Reconstruction from the Fourier Transform on surface with $R_0$ = {self.major_radius:.3f} at [t] = {t}")
        ax.set_xlabel(r"$\phi$ (Toroidal Angle)")
        ax.set_ylabel(r"$\theta$ (Poloidal Angle)")
        cf = ax.contourf(phi_grid,theta_grid,np.flip(B_n_ifft.real,axis=1),levels=50,vmin=B_n_ifft.real.min(),vmax=B_n_ifft.real.max(),cmap="viridis")
        fig.colorbar(cf,label="Minor Radial Magnetic Field")
        if show_plot:
            plt.show()
    
    def sort_fft_indices_and_mesh(self,B_n_fft,n_modes,m_modes):
        '''! Sort the Fourier Transformed B mesh and the poloidal and toroidal mode numbers generated by the FFT functions in this class while considering helicity

        @param B_n_fft The Fast Fourier Transformed B mesh, the first output of self.fft(B) [ntheta,nphi]
        @param n_modes The toroidal mode numbers, the second output of self.fft(B) [ntheta]
        @param m_modes The poloidal mode numbers, the third output of self.fft(B) [nphi]
        @result `B_n_sorted` The sorted FFT'ed B mesh [ntheta,nphi], `n_modes_sorted` The sorted toroidal mode numbers [ntheta], `m_modes_sorted` The sorted poloidal mode numbers [nphi]
        '''
    
        m_sorted_indices = np.argsort(m_modes)
        n_sorted_indices = np.argsort(n_modes)
        m_modes_sorted = m_modes[m_sorted_indices]
        n_modes_sorted = n_modes[n_sorted_indices]
        B_n_sorted = B_n_fft[:,n_sorted_indices][m_sorted_indices,:]

        return B_n_sorted, n_modes_sorted, m_modes_sorted
    
    def save_surfmn(self,t,filename,eliminate_negative_n=True,scale=1e-4,d_phi=None):
        '''! Sort the Fast Fourier Transformed B mesh in PEST or Hamada coordinate at time step `t` and write it to a SURFMN-formatted file

        @param t The time step at which the magnetic values are desired
        @param filename Filename of the file
        @param eliminate_negative_n Whether values of negative n modes should be written to the file
        @param scale The scaling of the values in the file
        @param d_phi Hamada phase shifts [ntheta]
        '''
            
        if d_phi is not None and d_phi[0] == d_phi[-1]:
            d_phi = d_phi[0:-1]

        B = np.roll(self.get_B_mesh(t)[:,::-1],shift=1,axis=1) if self.helicity == 1 else self.get_B_mesh(t)
        if d_phi is None:
            B_n_fft, n_modes, m_modes = self.fft2(B)
        elif len(d_phi) == len(self.radial_positions):
            B_n_fft, n_modes, m_modes = self.fft2(B,theta_coord='Hamada',d_phi=d_phi)
        else:
            raise ValueError('The d_phi input should have the same dimension as the radial/axial positions')
        
        B_n_sorted, n_modes_sorted, m_modes_sorted = self.sort_fft_indices_and_mesh(B_n_fft,n_modes,m_modes)

        # Prepare data
        mmax = m_modes_sorted[-1]
        mmin = m_modes_sorted[0]
        nmax = n_modes_sorted[-1]
        nmin = 0 if eliminate_negative_n else n_modes_sorted[0]
        mask = (n_modes_sorted>=nmin)
        B_fft = B_n_sorted[:,mask]

        num_m = mmax-mmin+1
        num_n = nmax-nmin+1
        dcosmn = np.zeros((num_m,num_n))
        dsinmn = np.zeros((num_m,num_n))

        print('Parameters to be set in GPEC:')
        print(f'nmin = {nmin}')
        print(f'nmax = {nmax}')
        print(f'mmin = {mmin}')
        print(f'mmax = {mmax}')

        # Write file
        for i in range(num_m):
            for j in range(num_n):
                val = B_fft[i,j]
                dcosmn[i,j] = val.real
                dsinmn[i,j] = val.imag

        with open(f"{filename}.dat","w") as f:
            if self.helicity == -1:
                mmin_h = mmax
                mmax_h = mmin-1
                step = -1
            else:
                mmin_h = mmin
                mmax_h = mmax+1
                step = 1
            for i in range(mmin_h,mmax_h,step):
                row_i = i-mmin

                # Write line for dcosmn(i, j=nmin...nmax)
                l = 0
                for j in range(nmin,nmax+1):
                    col_j = j-nmin
                    val_cos = dcosmn[row_i,col_j]/scale 
                    f.write(f"{val_cos:12.6f}")
                    l = l+1
                    if l == 25:
                        f.write("\n")
                        l = 0
                f.write("\n")

                # Write line for dsinmn(i, j=nmin...nmax)
                l = 0
                for j in range(nmin,nmax+1):
                    col_j = j-nmin
                    val_sin = dsinmn[row_i, col_j]/scale
                    f.write(f"{val_sin:12.6f}")
                    l = l+1
                    if l == 25:
                        f.write("\n")
                        l = 0
                f.write("\n")
    
    def plot_sensor_output(self,t,fig=None,ax=None):
        '''! Plot the magnetic field contour of the sensors on the surface at time step `t`

        @param t The time step during the time evolution
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        '''

        if hasattr(self,'hist_file') == False:
            raise AttributeError("The hist file of magnetic sensor values during time evolution has to be provided using the function load_histfile")
        else:
            phi_grid, theta_grid = np.meshgrid(np.linspace(0,2*np.pi,self.nphi,endpoint=False),self.theta_list)
            B_n = self.get_B_mesh(t)

            if ax is None and fig is None:
                fig,ax = plt.subplots(figsize=(10,6),constrained_layout=True)
                show_plot = True
            elif ax is not None and fig is not None:
                show_plot = False
            else:
                raise ValueError('For customized plotting, please provide both fig and ax as inputs.')
            ax.set_title(rf"Magnetic Field on surfaces of toroidal planes with $R_0$ = {self.major_radius:.3f} at [t] = {t}")
            ax.set_xlabel(r"$\phi$ (radians)")
            ax.set_ylabel(r"$\theta$ (radians)")
            cf = ax.contourf(phi_grid,theta_grid,np.flip(B_n,axis=1),vmax=B_n.max(),vmin=B_n.min(),levels=50,cmap="viridis")
            fig.colorbar(cf,label="Radial Magnetic Field (Tesla)")
            if show_plot:
                plt.show()

    def plot_sensor_output_on_surface(self,t,figs=None,axes=None):
        '''! Plot the normal magnetic field on the surface at time step `t`

        @param t The time step during the time evolution
        @param figa Maximumly two Matplotlib figures for plotting
        @param axes Two Matplotlib axes for plotting
        '''
    
        if axes is None and figs is None:
            fig1,ax1 = plt.subplots(constrained_layout=True)
            fig2,ax2 = plt.subplots(constrained_layout=True)
            figs = [fig1,fig2]
            axes = [ax1,ax2]
            show_plot = True
        elif axes is not None and figs is not None:
            axes = np.array(axes).flatten()
            figs = np.array(figs).flatten()
            if len(axes)!=2:
                raise ValueError('For customized plotting, two axes are required for plotting the surface at two cross sections.')
            if len(figs) == 1:
                figs = [figs[0],figs[0]]
            elif len(figs) != 2:
                raise ValueError('For customized plotting, maximumlly two figures are required.')
            show_plot = False
        else:
            raise ValueError('For customized plotting, please provide both figures and axes as inputs.')
        
        for i, angle in zip([0,self.nphi//2],np.linspace(2*np.pi,0,self.nphi,endpoint=False)[[0,self.nphi//2]]):
            points = np.array([self.radial_positions,self.axial_positions]).T.reshape(-1,1,2)
            points_closed = np.concatenate([points,points[:1]],axis=0)
            segments = np.concatenate([points_closed[:-1],points_closed[1:]],axis=1)
            z = self.get_B_mesh(t)[:,i]
            lc = LineCollection(segments,cmap='RdBu_r' if z.max()>=0 else 'RdBu',norm=plt.Normalize(z.min(),z.max()))
            lc.set_array(z[:-1])
            lc.set_linewidth(10)

            axes[i//(self.nphi//2)].add_collection(lc)
            axes[i//(self.nphi//2)].autoscale()
            figs[i//(self.nphi//2)].colorbar(lc,ax=axes[i//(self.nphi//2)],label='Radial Magnetic Field (Tesla)')
            axes[i//(self.nphi//2)].plot(np.concatenate([self.radial_positions,self.radial_positions[:1]]),np.concatenate([self.axial_positions,self.axial_positions[:1]]),color='black',linewidth=0.8,zorder=2)
            axes[i//(self.nphi//2)].set_xlabel('R')
            axes[i//(self.nphi//2)].set_ylabel('Z')
            axes[i//(self.nphi//2)].set_title(rf"Magnetic Field @ $\phi$ = {angle:.4f} at [t] = {t}")
            axes[i//(self.nphi//2)].set_aspect('equal')     
        if show_plot:
            plt.show()

    def plot_1D_fourier_amplitude(self,t,harmonics,axis,ax=None,d_phi=None,part='r'):
        '''! Plot real fourier amplitude of 1D Fast Fourier Transform of the magnetic mesh in `axis` at time step `t`

        @param t The time step during the time evolution
        @param harmonics Array of mode number in the `axis` dimension, whose real amplitudes are to be plotted against the other dimension [:]
        @param axis The direction in which the harmonics is to be obtained through 1D FFT; 1 is toroidal, and 0 is poloidal.
        @param ax Matplotlib axis for plotting
        @param d_phi Hamada phase shifts [ntheta]
        @param part Plot real ('r'), imaginary ('i'), or absolute ('a') amplitudes
        '''

        if axis == 1:
            B_n = self.get_B_mesh(t) if self.helicity == -1 else np.roll(self.get_B_mesh(t)[:,::-1],shift=1,axis=1)
            n_modes = np.fft.fftfreq(self.nphi)*self.nphi
            toroidal_harmonics=np.fft.fft(B_n,axis=1,norm="forward")
            # Apply phase shift
            if d_phi is not None and len(d_phi) == len(self.radial_positions):
                toroidal_harmonics *= np.exp(-1j*np.outer(d_phi,n_modes))
            elif d_phi is not None:
                raise ValueError('The d_phi input should have the same dimension as the radial/axial positions')
            toroidal_harmonics *= 2.0
            toroidal_harmonics[:,0] /= 2.0
            if len(n_modes)%2==0:
                toroidal_harmonics[:][len(n_modes)//2] /= 2.0
            if ax is None:
                _,ax = plt.subplots(figsize=(8,6),constrained_layout=True)
                show_plot = True
            else:
                show_plot = False
            for harmonic in harmonics:
                n_indice = np.where((n_modes == harmonic))[0]
                toroidal_mode = toroidal_harmonics[:,n_indice]
                if part =='r':
                    amplitude = toroidal_mode.real.flatten()
                    lb = 'Real'
                elif part == 'i':
                    amplitude = toroidal_mode.imag.flatten()
                    lb = 'Imag'
                elif part == 'a':
                    amplitude = abs(toroidal_mode).flatten()
                    lb = 'Abs'
                else:
                    raise ValueError("Input of 'part' is invalid.")
                ax.plot(self.theta_list,amplitude,label=f'n={harmonic}, {lb}')
            ax.legend()
            ax.set_title(f"1D FFT Amplitude of Toroidal Mode at [t] = {t}")
            ax.set_xlabel(r"$\theta$ (radians)")
            ax.set_ylabel(f"Mode amplitude (Tesla)")
            ax.grid()
            if show_plot:
                plt.show()
        elif axis == 0:
            B_n = self.get_B_mesh(t)
            m_modes = np.fft.fftfreq(self.ntheta)*self.ntheta
            poloidal_harmonics=np.fft.fft(B_n,axis=0,norm="forward")
            poloidal_harmonics *= 2.0
            poloidal_harmonics[0,:] /= 2.0
            if len(m_modes)%2==0:
                toroidal_harmonics[len(m_modes)//2][:] /= 2.0
            if ax is None:
                _,ax = plt.subplots(figsize=(8,6),constrained_layout=True)
                show_plot = True
            else:
                show_plot = False
            phi_list = np.linspace(0,2*np.pi,self.nphi,endpoint=False)
            for harmonic in harmonics:
                m_indice = np.where((m_modes == harmonic))[0]
                poloidal_mode = np.roll(poloidal_harmonics[m_indice,:].flatten()[::-1],shift=1)
                if part =='r':
                    amplitude = poloidal_mode.real.flatten()
                    lb = 'Real'
                elif part == 'i':
                    amplitude = poloidal_mode.imag.flatten()
                    lb = 'Imag'
                elif part == 'a':
                    amplitude = abs(poloidal_mode).flatten()
                    lb = 'Abs'
                else:
                    raise ValueError("Input of 'part' is invalid.")
                ax.plot(phi_list,poloidal_mode,label=f'm={harmonic}, {lb}')
            ax.legend()
            ax.set_title(f"1D FFT Amplitude of Poloidal Mode at [t] = {t}")
            ax.set_xlabel(r"$\phi$ (radians)")
            ax.set_ylabel(f"Mode amplitude (Tesla)")
            ax.grid()
            if show_plot:
                plt.show()
        else:
            raise ValueError('Unspecified or incorrect axis number, should be either 1 (toroidal harmonics) or 0 (poloidal harmonics)')
        
    def plot_m_over_n_amplitude(self,m_list,n_list,t_max,dt,t_min=0,ax=None,d_phi=None):
        '''! Plot the 2D Fast Fourier Transformed amplitude of the magnetic field for `m/n` (poloidal harmonic / toroidal harmonic) mode with respect to time

        @param m_list The list of poloidal modes to be visualize [:]
        @param n_list The list of toroidal modes to be visualize [:]
        @param t_max The ending time step of the plot
        @param dt The duration of each time step in unit of seconds
        @param t_min The starting time step of the plot
        @param ax Matplotlib axis for plotting
        @param d_phi Hamada phase shifts [ntheta]
        '''

        if (len(m_list)!=len(n_list)):
            raise ValueError('Input lists of poloidal and toroidal modes should have the same length')
        t_array = np.array(range(t_min,t_max+1))
        mode_amplitudes = np.zeros((len(m_list),len(t_array)))
        for t in t_array:
            B = self.get_B_mesh(t) if self.helicity == -1 else np.roll(self.get_B_mesh(t)[:,::-1],shift=1,axis=1)
            if d_phi is None:
                B_n_fft, n_modes, m_modes = self.fft2(B)
            elif len(d_phi) == len(self.radial_positions):
                B_n_fft, n_modes, m_modes = self.fft2(B,theta_coord='Hamada',d_phi=d_phi)
            else:
                raise ValueError('The d_phi input should have the same dimension as the radial/axial positions')
            B_n_sorted, n_modes_sorted, m_modes_sorted = self.sort_fft_indices_and_mesh(B_n_fft,n_modes,m_modes)

            i = 0
            for m,n in zip(m_list,n_list):
                m_indices = np.where((m_modes_sorted == m))[0]
                n_indices = np.where((n_modes_sorted == n))[0]
                mode_amplitudes[i][t] = B_n_sorted[m_indices,n_indices].real
                i+=1
    
        if ax is None:
            _,ax = plt.subplots(figsize=(8,6),constrained_layout=True)
            show_plot = True
        else:
            show_plot = False
        for j in range(mode_amplitudes.shape[0]):
            plt.plot(t_array*dt*1e3,mode_amplitudes[j],label=f'{m_list[j]}/{n_list[j]}')
        ax.set_ylabel('Mode amplitudes (Tesla)')
        ax.set_xlabel('Time (ms)')
        ax.set_title(r'Amplitude of $m/n$ modes in time')
        ax.legend()
        if show_plot:
            plt.show()

    def field_fourier_amplitude_contour(self,t,m_min,m_max,n_min,n_max,fig=None,ax=None,d_phi=None):
        '''! Plot the fourier amplitude of the magnetic field in the fourier space

        @param t The time step during the time evolution
        @param m_min The lower limit of the poloidal harmonics to be included in the contour
        @param m_max The upper limit of the poloidal harmonics to be included in the contour
        @param n_min The lower limit of the toroidal harmonics to be included in the contour
        @param n_max The upper limit of the toroidal harmonics to be included in the contour
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        @param d_phi Hamada phase shifts [ntheta]
        '''

        B = self.get_B_mesh(t) if self.helicity == -1 else np.roll(self.get_B_mesh(t)[:,::-1],shift=1,axis=1)
        if d_phi is None:
            B_n_fft, n_modes, m_modes = self.fft2(B)
        elif len(d_phi) == len(self.radial_positions):
            B_n_fft, n_modes, m_modes = self.fft2(B,theta_coord='Hamada',d_phi=d_phi)
        else:
            raise ValueError('The d_phi input should have the same dimension as the radial/axial positions')
        
        B_n_sorted, n_modes_sorted, m_modes_sorted = self.sort_fft_indices_and_mesh(B_n_fft,n_modes,m_modes)

        n_indices = np.where((n_modes_sorted == n_min) | (n_modes_sorted == n_max))[0]
        m_indices = np.where((m_modes_sorted == m_min) | (m_modes_sorted == m_max))[0]

        n_grid_fft, m_grid_fft = np.meshgrid(n_modes_sorted[n_indices[0]:n_indices[1]+1],m_modes_sorted[m_indices[0]:m_indices[1]+1])
        
        if ax is None and fig is None:
            fig,ax = plt.subplots(figsize=(10,6),constrained_layout=True)
            show_plot = True
        elif ax is not None and fig is not None:
            show_plot = False
        else:
            raise ValueError('For customized plotting, please provide both fig and ax as inputs.')
        ax.set_title(f"Real Fourier Amplitude of the Normal Magnetic Field at the Surface at [t] = {t}")
        ax.set_ylabel(r"$n$")
        ax.set_xlabel(r"$m$")
        cf = ax.contourf(m_grid_fft,n_grid_fft,B_n_sorted[m_indices[0]:m_indices[1]+1,n_indices[0]:n_indices[1]+1].real,levels=50,cmap="viridis")
        fig.colorbar(cf,label=f"Mode real amplitude (Tesla)")
        if show_plot:
            plt.show()

    def plot_2D_fourier_amplitude(self,t,harmonics,toroidal_harmonics=True,d_phi=None,ax=None,x_type='modes',x_mode_min=None,x_mode_max=None):
        '''! Plot the 2D Fast Fourier Transformed amplitude of the mesh of magnetic values against poloidal/toroidal harmonics/angles

        @param t The time step during the time evolution
        @param harmonics List of (poloidal/toroidal) modes whose amplitudes are to be visualized in y axis [:]
        @param toroidal_harmonics Whether the input `harmonics` is toroidal or poloidal harmonics
        @param d_phi Hamada phase shifts [ntheta]
        @param ax Matplotlib axis for plotting
        @param x_type The variable on x axis ('modes' or 'angles')
        @param x_mode_min The min of the (toroidal/poloidal) mode number to be visualized on x axis
        @param x_mode_max The max of the (toroidal/poloidal) mode number to be visualized on x axis
        '''

        harmonics = np.array([harmonics]).flatten()
        if d_phi is not None and d_phi[0] == d_phi[-1]:
            d_phi = d_phi[0:-1]
        if x_type not in ['modes','angles']:
            raise ValueError("Unsupported x variable is provided. Accepts 'modes' and 'angles' only.")
        elif x_type == 'modes' and (x_mode_min is None or x_mode_max is None):
            raise ValueError('For x_type == "modes", both mode_min and mode_max should be provided.')
            
        B = np.roll(self.get_B_mesh(t)[:,::-1],shift=1,axis=1) if self.helicity == 1 else self.get_B_mesh(t)
        if d_phi is None:
            B_n_fft, n_modes, m_modes = self.fft2(B)
        elif len(d_phi) == len(self.radial_positions):
            B_n_fft, n_modes, m_modes = self.fft2(B,theta_coord='Hamada',d_phi=d_phi)
        else:
            raise ValueError('The d_phi input should have the same dimension as the radial/axial positions')
        
        if x_type == 'modes':
            B_n_sorted, n_modes_sorted, m_modes_sorted = self.sort_fft_indices_and_mesh(B_n_fft,n_modes,m_modes)
            if ax is None:
                _,ax = plt.subplots(figsize=(20,10),dpi=200,constrained_layout=True)
                show_plot = True
            else:
                show_plot = False
            cmap = plt.get_cmap('tab10')
            if toroidal_harmonics:
                m_range = np.where((m_modes_sorted == x_mode_min) | (m_modes_sorted == x_mode_max))[0]
                mode_idx = [np.where(n_modes_sorted == harmonic)[0][0] for harmonic in harmonics]
                
                for i, idx in enumerate(mode_idx):
                    color = cmap(i/len(mode_idx))
                    ax.plot(m_modes_sorted[m_range[0]:m_range[1]+1],B_n_sorted[m_range[0]:m_range[1]+1,idx].real,color=color,label=f"n={harmonics[i]}, real")
                    ax.plot(m_modes_sorted[m_range[0]:m_range[1]+1],B_n_sorted[m_range[0]:m_range[1]+1,idx].imag,linestyle='--',color=color,label=f"n={harmonics[i]}, imag")
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Toroidal Modes at [t] = {t}",fontsize=20)
                ax.set_xlabel(f"Poloidal Harmonics ($m$)")
                ax.set_ylabel(f"Mode Amplitudes (Tesla)")
                ax.set_xticks(range(x_mode_min,x_mode_max+1))
                ax.grid()
                if show_plot:
                    plt.show()
            else:
                n_range = np.where((n_modes_sorted == x_mode_min) | (n_modes_sorted == x_mode_max))[0]
                mode_idx = [np.where(m_modes_sorted == harmonic)[0][0] for harmonic in harmonics]

                for i, idx in enumerate(mode_idx):
                    color = cmap(i/len(mode_idx))
                    ax.plot(n_modes_sorted[n_range[0]:n_range[1]+1],B_n_sorted[idx,n_range[0]:n_range[1]+1].real,color=color,label=f"m={harmonics[i]}, real")
                    ax.plot(n_modes_sorted[n_range[0]:n_range[1]+1],B_n_sorted[idx,n_range[0]:n_range[1]+1].imag,linestyle='--',color=color,label=f"m={harmonics[i]}, imag")
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Poloidal Modes at [t] = {t}",fontsize=20)
                ax.set_xlabel(f"Toroidal Harmonics ($n$)")
                ax.set_ylabel(f"Mode Amplitudes (Tesla)")
                ax.set_xticks(range(x_mode_min,x_mode_max+1))
                ax.grid()
                if show_plot:
                    plt.show()
        else:
            if ax is None:
                _,ax = plt.subplots(figsize=(10,6),constrained_layout=True)
                show_plot = True
            else:
                show_plot = False
            if toroidal_harmonics:
                mode_idx = [np.where(n_modes == harmonic)[0][0] for harmonic in harmonics]
                for i, idx in enumerate(mode_idx):
                    ax.plot(self.theta_list,B_n_fft[:,idx].real,label=f"n={harmonics[i]}, real")
                    ax.plot(self.theta_list,B_n_fft[:,idx].imag,linestyle='--',label=f"n={harmonics[i]}, imag")
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Toroidal Modes at [t] = {t}")
                ax.set_xlabel(r"$\theta$ (radians)")
                ax.set_ylabel(f"Mode Amplitudes (Tesla)")
                if show_plot:
                    plt.show()
            else:
                mode_idx = [np.where(m_modes == harmonic)[0][0] for harmonic in harmonics]
                phi_list = np.linspace(0,2*np.pi,self.nphi,endpoint=False)
                for i, idx in enumerate(mode_idx):
                    if self.helicity == -1:
                        ax.plot(phi_list,np.roll(B_n_fft[idx,:].real[::-1],shift=1),label=f"m={harmonics[i]}, real")
                        ax.plot(phi_list,np.roll(B_n_fft[idx,:].imag[::-1],shift=1),linestyle='--',label=f"m={harmonics[i]}, imag")
                    else:
                        ax.plot(phi_list,B_n_fft[idx,:].real,label=f"m={harmonics[i]}, real")
                        ax.plot(phi_list,B_n_fft[idx,:].imag,linestyle='--',label=f"m={harmonics[i]}, imag")
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Poloidal Modes at [t] = {t}")
                ax.set_xlabel(r"$\phi$ (radians)")
                ax.set_ylabel(f"Mode Amplitudes (Tesla)")
                if show_plot:
                    plt.show()

    def plot_sensor_signal_against_angle(self,t,theta=True,ax=None):
        '''! Plot the value of normal magnetic fields over theta/phi at at phi/theta = 0

        @param t The time step during the time evolution
        @param theta Plot against theta (True) or phi (False)
        @param ax Matplotlib axis for plotting
        '''

        B_n = self.get_B_mesh(t)
        if ax is None:
            _,ax = plt.subplots(figsize=(10,6),constrained_layout=True)
            show_plot = True
        else:
            show_plot = False
        if theta:
            ax.plot(self.theta_list/2/np.pi, B_n[:,0])
            ax.set_title(rf"Magnetic Field on surface @ $\phi$=0 at [t] = {t}")
            ax.set_xlabel(r"$\theta$ (radians)")
            ax.set_ylabel(f"Magnetic Field (Tesla)")
            if show_plot:
                plt.show()
        else:
            phi_list = np.linspace(0,2*np.pi,self.nphi,endpoint=False)
            ax.plot(phi_list/2/np.pi, np.flip(B_n[0,:],axis=1))
            ax.set_title(rf"Magnetic Field on surface @ $\theta$=0 at [t] = {t}")
            ax.set_xlabel(r"$\phi$ (radians)")
            ax.set_ylabel(f"Magnetic Field (Tesla)")
            if show_plot:
                plt.show()
