'''! Python interface for ThinCurr thin-wall eddy current functionality

@authors Thomas Wang
@date June 2025
@ingroup doxy_oft_python
'''
import re
import xml.etree.ElementTree as ET
import numpy as np
from .sensor import Mirnov, save_sensors
from ..io import histfile
from ..util import mu0

class torus_fourier_sensor():
    '''! Interface that sets up normal magnetic field probes at a desired surface, Fourier analyze the signals, and output them for GPEC'''

    def __init__(self,radial_positions,axial_positions,major_radius,helicity,hamada_dphi=None):
        '''! Initialize the interface with surface

        @param radial_positions The R coordinates that defines the surface [ntheta]
        @param axial_positions The Z coordinates that defines the surface [ntheta]
        @param major_radius The major radius of the device
        @param helicity Sign of the plasma helicity
        @param hamada_dphi Toroidal angle minus magnetic angle at each surface point,
               in radians [ntheta]. Trimmed and sorted together with the surface points;
               methods that accept a `hamada_dphi` argument use this value when not
               given one.
        '''

        if helicity not in [1,-1]:
            raise ValueError('Helicity should be either 1 or -1')
        self.major_radius = major_radius
        self.helicity = helicity
        npts_in = len(radial_positions)
        if (abs(radial_positions[0] - radial_positions[-1]) < 1e-13) or (abs(axial_positions[0] - axial_positions[-1]) < 1e-13):
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
        self.ntheta = len(self.radial_positions)
        if hamada_dphi is None:
            self.hamada_dphi = None
        else:
            hamada_dphi = np.asarray(hamada_dphi)
            if len(hamada_dphi) == npts_in and npts_in != self.ntheta:
                hamada_dphi = hamada_dphi[0:-1]   # drop the duplicate endpoint, as for R,Z
            if len(hamada_dphi) != self.ntheta:
                raise ValueError('`hamada_dphi` should have the same dimension as the '
                                 'radial/axial positions.')
            # the surface is reordered by poloidal angle above; keep dphi aligned with it
            self.hamada_dphi = hamada_dphi[sorted_indices]
        self.hist_file = None

    @classmethod
    def from_eqdsk(cls,radial_positions,axial_positions,eqdsk_file,center='raxis',hamada_dphi=None,verbose=True):
        '''! Construct the interface, taking `major_radius` and `helicity` from an equilibrium

        Avoids hand-entering values that are already defined by the equilibrium. The
        g-file is normalized to COCOS 1 on read, so the derived helicity does not depend
        on the sign conventions used to write the file.

        @param radial_positions The R coordinates that defines the surface [ntheta]
        @param axial_positions The Z coordinates that defines the surface [ntheta]
        @param eqdsk_file Path to a gEQDSK/eqdsk equilibrium file
        @param center Origin for the poloidal angle: 'raxis' (magnetic axis, default),
               'rcentr' (vacuum field reference radius), or 'geometric'
               (midpoint of the supplied surface, uses no equilibrium geometry)
        @param hamada_dphi Toroidal angle minus magnetic angle at each surface point [ntheta]
        @param verbose Print the derived values so they stay visible
        @result A `torus_fourier_sensor` instance
        '''
        # imported here rather than at module scope to keep ThinCurr independent of TokaMaker
        from ..TokaMaker.eqdsk import read_geqdsk
        eq = read_geqdsk(eqdsk_file)

        if center == 'raxis':
            major_radius = float(eq.R_mag)
        elif center == 'rcentr':
            major_radius = float(eq.R_center)
        elif center == 'geometric':
            major_radius = 0.5*(np.max(radial_positions)+np.min(radial_positions))
        else:
            raise ValueError("center should be one of 'raxis', 'rcentr' or 'geometric'")

        # Helicity is the sign of the product of toroidal field and plasma current. The
        # product is invariant under COCOS conversion and under reversing both Bt and Ip,
        # so a field-reversed AND current-reversed equilibrium is still helicity +1.
        b_center = float(eq.B_center)
        plasma_current = float(eq.Ip)
        signed = b_center*plasma_current
        if signed == 0.0:
            raise ValueError('Cannot determine helicity: B_center*Ip is zero in %s' % eqdsk_file)
        helicity = 1 if signed > 0.0 else -1

        if verbose:
            print("torus_fourier_sensor.from_eqdsk(%s):" % eqdsk_file)
            print("    major_radius = %.7f  (%s)" % (major_radius,center))
            print("    helicity     = %+d       (B_center=%+.4f, Ip=%+.4g, COCOS %d)"
                  % (helicity,b_center,plasma_current,eq.cocos))

        return cls(radial_positions,axial_positions,major_radius,helicity,hamada_dphi=hamada_dphi)

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
        @result line The line object of the surface with outward normal vectors plotted on `ax` if it is provided
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
        self.nphi = nphi
        for i, phi in enumerate(np.linspace(2*np.pi,0,nphi,endpoint=False)):
            for j in range(self.ntheta):
                sensors.append(Mirnov(np.r_[self.radial_positions[j]*np.cos(phi),self.radial_positions[j]*np.sin(phi),self.axial_positions[j]],np.r_[outward_normals[j,0]*np.cos(phi),outward_normals[j,0]*np.sin(phi),outward_normals[j,1]],'B_{0}_{1}'.format(i+1,j+1)))

        save_sensors(sensors,filename=filename)

        if ax is not None:
            if len(np.array(ax).flatten())!=1:
                raise ValueError('For customized plotting, please provide a single axis for plotting the surface with outward normal vectors.')
            line, = ax.plot(np.concatenate([self.radial_positions,self.radial_positions[:1]]),np.concatenate([self.axial_positions,self.axial_positions[:1]]))
            ax.quiver(self.radial_positions,self.axial_positions,outward_normals[:,0],outward_normals[:,1],color='red',scale=20)
            ax.set_xlabel("R (m)")
            ax.set_ylabel("Z (m)")
            ax.set_ylim((-1.25*max(self.axial_positions),1.25*max(self.axial_positions)))
            ax.set_xlim((min(self.radial_positions)-0.1,max(self.radial_positions)+0.1))
            ax.set_title("Surface with Outward Normal Vectors")
            ax.grid(True)
            ax.set_aspect('equal')
            return line

    def load_histfile(self, hist_file_path='floops.hist', norm=None):
        '''! Loading histfile containing magnetic values at the surface collected by sensors

            @param hist_file_path Path to the histfile containing sensor values
            @param norm Normalization factor to be applied to the sensor values
        '''

        hist_file = histfile(hist_file_path)
        if hist_file.nfields-1 != self.ntheta*self.nphi:
            raise ValueError("The hist file provided might not be the output of the sensors defined in the current instance of the sensor_interface class.")
        else:
            if norm is not None:
                for name in hist_file._data:
                    if name != 'time':
                        hist_file._data[name] /= norm
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
    
    def get_complex_sensor_mesh_freq_response(self, probe_signals=None, filepath=None):
        '''! Generates the complex sensor mesh of shape (ntheta, nphi)
        
        @param probe_signals The array calculated from the result of `compute_freq_response` [2, ntheta * nphi]
        @param filepath Path to the text file (if reading from saved output)
        @result `complex_sensor_mesh` Numpy array of the complex amplitudes for frequency response [ntheta, nphi]
        '''
        
        # --- PATH 1: Using the direct array from the simulation ---
        if probe_signals is not None:
            # Combine Real (row 0) and Imaginary (row 1) into a 1D complex array
            complex_signals = probe_signals[0, :] + 1j * probe_signals[1, :]
            
        # --- PATH 2: Reading from the text file ---
        elif filepath is not None:
            with open(filepath, 'r') as f:
                lines = f.readlines()
                
            complex_signals = []
            for line in lines:
                if "Real" in line and "Imaginary" in line:
                    parts = line.replace(',', '').split()
                    real_val = float(parts[1])
                    imag_val = float(parts[3])
                    complex_signals.append(complex(real_val, imag_val))
                    
            complex_signals = np.array(complex_signals)
            
        else:
            raise ValueError("You must provide either 'probe_signals' (from frequency response calculation) or the 'filepath' to its storage.")

        # --- SAFETY CHECK ---
        expected_size = self.nphi * self.ntheta
        if len(complex_signals) != expected_size:
            raise ValueError(f"Data size ({len(complex_signals)}) does not match expected nphi*ntheta ({expected_size})")

        # --- RESHAPE AND RETURN ---
        # The data is ordered by phi first, then theta. 
        # We reshape to (nphi, ntheta) and transpose to (ntheta, nphi)
        complex_sensor_mesh = complex_signals.reshape((self.nphi, self.ntheta)).T
        
        return complex_sensor_mesh
        
    def fft2(self,B,hamada_dphi=None):
        '''! Performs a 2D Fast Fourier Transform in PEST/Hamada coordinates on the magnetic field matrix probed by sensors with proper Nyquist handling

        @param B Input matrix of shape [ntheta,nphi]
        @param hamada_dphi Hamada phase shifts [ntheta]
        @result `B_n` The Fast Fourier Transformed matrix [ntheta,nphi], `n_modes` The toroidal mode numbers [nphi], `m_modes` The poloidal mode numbers [ntheta]
        '''
        if hamada_dphi is None:
            hamada_dphi = self.hamada_dphi

        B = np.roll(B[:,::-1],shift=1,axis=1) if self.helicity == 1 else B

        n_theta, n_phi = B.shape
        n_modes = np.round(np.fft.fftfreq(n_phi)*n_phi).astype(int)
        m_modes = np.round(np.fft.fftfreq(n_theta)*n_theta).astype(int)

        if hamada_dphi is None:
            B_n = np.fft.fft2(B,norm='forward')
        elif hamada_dphi is not None and len(hamada_dphi) == len(self.radial_positions):
            B_n_toroidal = np.fft.fft(B,axis=1,norm='forward')
            phase_matrix = np.exp(-1j*np.outer(hamada_dphi,n_modes))
            B_n_toroidal_shifted = B_n_toroidal * phase_matrix
            B_n = np.fft.fft(B_n_toroidal_shifted,axis=0,norm='forward')
        else:
            raise ValueError('`hamada_dphi` should have the same dimension as the radial/axial positions.')

        B_n *= 2.0
        B_n[0,0] /= 2.0
        if n_theta%2==0:
            B_n[n_theta//2,0] /= 2.0
        if n_phi%2==0:
            B_n[0,n_phi//2] /= 2.0
        if n_theta%2==0 and n_phi%2==0:
            B_n[n_theta//2,n_phi//2] /=2.0

        return B_n, n_modes, m_modes

    def ifft2(self,B_n,hamada_dphi=None):
        '''! Performs an inverse 2D Fast Fourier Transform in PEST/Hamada coordinates on the transformed magnetic field matrix with proper Nyquist handling

        @param B_n Input FFT'ed matrix of shape [ntheta,nphi]
        @param hamada_dphi Hamada phase shifts [ntheta]
        @result `B_ifft` The reconstructed matrix [ntheta,nphi]
        '''
        if hamada_dphi is None:
            hamada_dphi = self.hamada_dphi

        n_theta, n_phi = B_n.shape

        B_n /= 2.0
        B_n[0,0] *= 2.0
        if n_theta%2 == 0:
            B_n[n_theta//2][0] *= 2.0
        if n_phi%2 == 0:
            B_n[0][n_phi//2] *= 2.0
        if n_theta%2 == 0 and n_phi%2 == 0:
            B_n[n_theta//2][n_phi//2] *=2.0

        if hamada_dphi is None:
            return np.fft.ifft2(B_n,norm="forward")
        elif hamada_dphi is not None and len(hamada_dphi) == len(self.radial_positions):
            B_ifft_poloidal_shifted = np.fft.ifft(B_n,axis=0,norm='forward')
            n_modes = np.round(np.fft.fftfreq(n_phi)*n_phi).astype(int)
            inverse_phase = np.exp(1j*np.outer(hamada_dphi,n_modes))
            B_ifft_poloidal = B_ifft_poloidal_shifted*inverse_phase
            return np.fft.ifft(B_ifft_poloidal,axis=1,norm='forward')
        else:
            raise ValueError('`hamada_dphi` should have the same dimension as the radial/axial positions.')

    def plot_inverse_2D_fourier_transform(self,t,fig,ax):
        '''! Use ifft2 to reconstruct the sensor signal at time step `t`

        @param t The time step of the sensor signals at which fft and ifft is to be performed
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        @result cf The contour plot of the reconstructed sensor signals on the surface at time step `t`
        @result cbar The colorbar of the contour plot
        '''

        B_n_fft, _, _ = self.fft2(self.get_B_mesh(t))
        B_n_ifft = self.ifft2(B_n_fft)
        # print(f'The largest difference between the original and reconstructed sensor signals is {np.linalg.norm(abs(B_n_ifft.real-self.get_B_mesh(t).real),np.inf)}')
        phi_grid, theta_grid = np.meshgrid(np.linspace(0,2*np.pi,self.nphi,endpoint=False),self.theta_list)
        if len(np.array(ax).flatten())!=1 or len(np.array(fig).flatten())!=1:
            raise ValueError('For customized plotting, please provide a single figure and a single axis for plotting the a contour.')
        ax.set_title(fr"Reconstruction from the Fourier Transform on surface with $R_0$ = {self.major_radius:.3f} at [t] = {t}")
        ax.set_xlabel(r"$\phi$ (Toroidal Angle)")
        ax.set_ylabel(r"$\theta$ (Poloidal Angle)")
        cf = ax.contourf(phi_grid,theta_grid,np.flip(B_n_ifft.real,axis=1),levels=50,vmin=B_n_ifft.real.min(),vmax=B_n_ifft.real.max(),cmap="viridis")
        cbar = fig.colorbar(cf,label="Minor Radial Magnetic Field")
        cbar.ax.ticklabel_format(style='sci', scilimits=(-3, 3))
        return cf, cbar

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
    
    def save_spectrum(self,t,filename,eliminate_negative_n=True,scale=1e-4,hamada_dphi=None,sensor_mesh=None,data_type='surfmn'):
        '''! Sort the Fast Fourier Transformed B mesh in PEST or Hamada coordinate at time step `t` and write it to a SURFMN-formatted file

        @param t The time step at which the magnetic values are desired
        @param filename Filename of the file
        @param eliminate_negative_n Whether values of negative n modes should be written to the file
        @param scale The scaling of the values in the file (ignored for `data_type='vac3d'`, which GPEC reads with scale=1.0, i.e. Tesla)
        @param hamada_dphi Hamada phase shifts [ntheta]
        @param sensor_mesh Customized sensor signal mesh to be written to the file (often from frequency response calculation) [ntheta, nphi]
        @param data_type GPEC `data_type` the file is written for: 'surfmn' (fixed-point Gauss, `(1x,25f12.6)`, ~6 decimals absolute)
               or 'vac3d' (scientific Tesla, `(11(1x,e15.8))`, ~9 significant digits). GPEC parses both into the same
               control-surface spectrum; prefer 'vac3d' when file precision matters.
        '''
        if hamada_dphi is None:
            hamada_dphi = self.hamada_dphi

        if sensor_mesh is None:
            B = self.get_B_mesh(t)
        elif sensor_mesh.shape != (self.ntheta,self.nphi):
            raise ValueError("The shape of sensor_mesh does not match the sensor dimensions.")
        else:
            B = sensor_mesh

        if hamada_dphi is None:
            B_n_fft, n_modes, m_modes = self.fft2(B)
        else:
            if abs(hamada_dphi[0] - hamada_dphi[-1]) < 1e-10:
                hamada_dphi = hamada_dphi[0:-1]
            B_n_fft, n_modes, m_modes = self.fft2(B,hamada_dphi=hamada_dphi)

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
        print(f'data_type = "{data_type}"')
        print(f'nmin = {nmin}')
        print(f'nmax = {nmax}')
        print(f'mmin = {mmin}')
        print(f'mmax = {mmax}')

        if data_type == 'surfmn':
            fmt, per_line = "{:12.6f}", 25
        elif data_type == 'vac3d':
            fmt, per_line = " {:15.8E}", 11
            scale = 1.0  # GPEC's vac3d branch keeps scale=1.0 (fields in Tesla)
        else:
            raise ValueError("data_type must be 'surfmn' or 'vac3d'")

        # Write file
        for i in range(num_m):
            for j in range(num_n):
                val = B_fft[i,j]
                dcosmn[i,j] = val.real
                dsinmn[i,j] = val.imag

        with open(f"{filename}.dat","w") as f:
            # GPEC reverses its m read-order for helicity<0 only in the surfmn branch;
            # the vac3d branch always reads m ascending.
            if self.helicity == -1 and data_type == 'surfmn':
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
                count = 0
                for j in range(nmin,nmax+1):
                    col_j = j-nmin
                    val_cos = dcosmn[row_i,col_j]/scale
                    f.write(fmt.format(val_cos))
                    count = count+1
                    if count == per_line:
                        f.write("\n")
                        count = 0
                f.write("\n")

                # Write line for dsinmn(i, j=nmin...nmax)
                count = 0
                for j in range(nmin,nmax+1):
                    col_j = j-nmin
                    val_sin = dsinmn[row_i, col_j]/scale
                    f.write(fmt.format(val_sin))
                    count = count+1
                    if count == per_line:
                        f.write("\n")
                        count = 0
                f.write("\n")

    def plot_sensor_output(self,t,fig,ax):
        '''! Plot the magnetic field contour of the sensors on the surface at time step `t`

        @param t The time step during the time evolution
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        @result cf The contour plot of the sensor values on the surface at time step `t`
        @result cbar The colorbar of the contour plot
        '''

        if not hasattr(self,'hist_file'):
            raise AttributeError("The hist file of magnetic sensor values during time evolution has to be provided using the function load_histfile")
        else:
            phi_grid, theta_grid = np.meshgrid(np.linspace(0,2*np.pi,self.nphi,endpoint=False),self.theta_list)
            B_n = self.get_B_mesh(t)

            if len(np.array(ax).flatten())!=1 or len(np.array(fig).flatten())!=1:
                raise ValueError('For customized plotting, please provide a single figure and a single axis for plotting the a contour.')
            ax.set_title(rf"Magnetic Field on surfaces of toroidal planes with $R_0$ = {self.major_radius:.3f} at [t] = {t}")
            ax.set_xlabel(r"$\phi$ (radians)")
            ax.set_ylabel(r"$\theta$ (radians)")
            cf = ax.contourf(phi_grid,theta_grid,np.flip(B_n,axis=1),vmax=B_n.max(),vmin=B_n.min(),levels=50,cmap="RdBu_r")
            cbar = fig.colorbar(cf,label="Minor Radial Magnetic Field (Tesla)")
            cbar.ax.ticklabel_format(style='sci', scilimits=(-3, 3))
            return cf, cbar

    def plot_sensor_output_on_surface(self,t,figs,axes):
        '''! Plot the normal magnetic field on the surface at time step `t`

        @param t The time step during the time evolution
        @param figs Maximumly two Matplotlib figures for plotting
        @param axes Two Matplotlib axes for plotting
        @result lc_list The list of LineCollection objects for the two cross sections
        @result boundary_list The list of surface boundary lines
        @result cbar_list The list of colorbars for the LineCollection objects
        '''

        axes = np.array(axes).flatten()
        figs = np.array(figs).flatten()
        if len(axes)!=2:
            raise ValueError('For customized plotting, two axes are required for plotting the surface at two cross sections.')
        if len(figs) == 1:
            figs = [figs[0],figs[0]]
        elif len(figs) != 2:
            raise ValueError('For customized plotting, maximumlly two figures are required.')

        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        lc_list = []
        boundary_list = []
        cbar_list = []
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
            cbar = figs[i//(self.nphi//2)].colorbar(lc,ax=axes[i//(self.nphi//2)],label='Normal Magnetic Field (Tesla)')
            cbar.ax.ticklabel_format(style='sci', scilimits=(-3, 3))
            boundary = axes[i//(self.nphi//2)].plot(np.concatenate([self.radial_positions,self.radial_positions[:1]]),np.concatenate([self.axial_positions,self.axial_positions[:1]]),color='black',linewidth=0.8,zorder=2)
            axes[i//(self.nphi//2)].set_xlabel('R')
            axes[i//(self.nphi//2)].set_ylabel('Z')
            axes[i//(self.nphi//2)].set_title(rf"Magnetic Field @ $\phi$ = {angle:.4f} at [t] = {t}")
            axes[i//(self.nphi//2)].set_aspect('equal')
            lc_list.append(lc)
            boundary_list.append(boundary)
            cbar_list.append(cbar)
        return lc_list, boundary_list, cbar_list

    def plot_complex_signal_amplitude(self, complex_sensor_mesh, fig, ax):
        '''! Plot the amplitude contour of the complex sensor values

        @param complex_sensor_mesh The complex sensor mesh generated by `get_complex_sensor_mesh_freq_response` [ntheta, nphi]
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        @result cf The contour plot of the sensor values on the surface at time step `t`
        @result cbar The colorbar of the contour plot
        '''
        amplitude = np.abs(complex_sensor_mesh)
        phi_grid, theta_grid = np.meshgrid(np.linspace(0, 2*np.pi, self.nphi, endpoint=False), self.theta_list)

        ax.set_title(rf"Frequency Response - Magnetic Field Amplitude ($R_0$ = {self.major_radius:.3f} m)")
        ax.set_xlabel(r"$\phi$ (radians)")
        ax.set_ylabel(r"$\theta$ (radians)")
        
        cf = ax.contourf(phi_grid, theta_grid, np.flip(amplitude, axis=1), 
                         vmax=amplitude.max(), vmin=amplitude.min(), levels=50, cmap="viridis")
        
        cbar = fig.colorbar(cf, ax=ax, label="Minor Radial Magnetic Field (Tesla)")
        cbar.ax.ticklabel_format(style='sci', scilimits=(-3, 3))
        
        return cf, cbar

    def plot_complex_signal_phase(self, complex_sensor_mesh, fig, ax):
        '''! Plot the phase contour of the complex sensor values

        @param complex_sensor_mesh The complex sensor mesh generated by `get_complex_sensor_mesh_freq_response` [ntheta, nphi]
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        @result cf The contour plot of the sensor values on the surface at time step `t`
        @result cbar The colorbar of the contour plot
        '''
        # np.angle returns the phase shift in radians from -pi to pi
        phase = np.angle(complex_sensor_mesh)
        phi_grid, theta_grid = np.meshgrid(np.linspace(0, 2*np.pi, self.nphi, endpoint=False), self.theta_list)

        ax.set_title(rf"Frequency Response - Magnetic Field Phase Shift ($R_0$ = {self.major_radius:.3f} m)")
        ax.set_xlabel(r"$\phi$ (radians)")
        ax.set_ylabel(r"$\theta$ (radians)")
        
        # 'twilight' is a cyclic colormap, meaning -pi and +pi are the same color, 
        # which is physically accurate for angles!
        # cf = ax.contourf(phi_grid, theta_grid, np.flip(phase, axis=1), 
        #                  vmax=np.pi, vmin=-np.pi, levels=50, cmap="twilight")
        cf = ax.pcolormesh(phi_grid, theta_grid, np.flip(phase, axis=1), 
                   vmax=np.pi, vmin=-np.pi, cmap="twilight", shading='nearest')
        
        cbar = fig.colorbar(cf, ax=ax, label="Phase (radians)")
        
        return cf, cbar 

    def plot_1D_fourier_amplitude(self,t,harmonics,ax,toroidal_harmonics=True,hamada_dphi=None,part='r'):
        '''! Plot real fourier amplitude of 1D Fast Fourier Transform of the magnetic mesh in `axis` at time step `t`

        @param t The time step during the time evolution
        @param harmonics Array of mode number in the `axis` dimension, whose real amplitudes are to be plotted against the other dimension [:]
        @param ax Matplotlib axis for plotting
        @param toroidal_harmonics Whether the input `harmonics` is toroidal or poloidal harmonics
        @param hamada_dphi Hamada phase shifts [ntheta]
        @param part The part of the fourier amplitude to be plotted ('r' for real, 'i' for imaginary, 'a' for absolute)
        @result line_list The list of line objects for each harmonic in `harmonics`
        @result FFTed_mesh The 1D FFTed mesh in the direction desired [ntheta,nphi]
        '''
        if hamada_dphi is None:
            hamada_dphi = self.hamada_dphi

        harmonics = np.array([harmonics]).flatten()
        if toroidal_harmonics:
            B_n = self.get_B_mesh(t) if self.helicity == -1 else np.roll(self.get_B_mesh(t)[:,::-1],shift=1,axis=1)
            n_modes = np.fft.fftfreq(self.nphi)*self.nphi
            toroidal_harmonics=np.fft.fft(B_n,axis=1,norm="forward")
            # Apply phase shift
            if hamada_dphi is not None:
                if abs(hamada_dphi[0] - hamada_dphi[-1]) < 1e-10:
                    hamada_dphi = hamada_dphi[0:-1]
                if len(hamada_dphi) == len(self.radial_positions):
                    toroidal_harmonics *= np.exp(-1j*np.outer(hamada_dphi,n_modes))
                else:
                    raise ValueError('The hamada_dphi input should have the same dimension as the radial/axial positions')
            toroidal_harmonics *= 2.0
            toroidal_harmonics[:,0] /= 2.0
            if len(n_modes)%2==0:
                toroidal_harmonics[:][len(n_modes)//2] /= 2.0
            line_list = []
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
                line = ax.plot(self.theta_list,amplitude,label=f'n={harmonic}, {lb}')
                line_list.append(line)
            ax.legend()
            ax.set_title(f"1D FFT Amplitude of Toroidal Mode at [t] = {t}")
            ax.set_xlabel(r"$\theta$ (radians)")
            ax.set_ylabel("Mode amplitude (Tesla)")
            ax.grid()
            ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
            return line_list, toroidal_harmonics
        else:
            B_n = self.get_B_mesh(t)
            m_modes = np.fft.fftfreq(self.ntheta)*self.ntheta
            poloidal_harmonics=np.fft.fft(B_n,axis=0,norm="forward")
            poloidal_harmonics *= 2.0
            poloidal_harmonics[0,:] /= 2.0
            if len(m_modes)%2==0:
                poloidal_harmonics[len(m_modes)//2][:] /= 2.0
            poloidal_harmonics = np.roll(poloidal_harmonics[:,::-1],shift=1,axis=1)
            phi_list = np.linspace(0,2*np.pi,self.nphi,endpoint=False)
            line_list = []
            for harmonic in harmonics:
                m_indice = np.where((m_modes == harmonic))[0]
                poloidal_mode = poloidal_harmonics[m_indice,:]
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
                line = ax.plot(phi_list,poloidal_mode,label=f'm={harmonic}, {lb}')
                line_list.append(line)
            ax.legend()
            ax.set_title(f"1D FFT Amplitude of Poloidal Mode at [t] = {t}")
            ax.set_xlabel(r"$\phi$ (radians)")
            ax.set_ylabel("Mode amplitude (Tesla)")
            ax.grid()
            ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
            return line_list, poloidal_harmonics

    def plot_m_over_n_amplitude(self,m_list,n_list,t_max,dt,ax,t_min=0,hamada_dphi=None):
        '''! Plot the 2D Fast Fourier Transformed amplitude of the magnetic field for `m/n` (poloidal harmonic / toroidal harmonic) mode with respect to time

        @param m_list The list of m of the m/n modes to be visualized [:]
        @param n_list The (list of) n of the m/n modes to be visualized [:]
        @param t_max The ending time step of the plot
        @param dt The duration of each time step in unit of seconds
        @param ax Matplotlib axis for plotting
        @param t_min The starting time step of the plot
        @param hamada_dphi Hamada phase shifts [ntheta]
        @result line_list The list of line objects for each m/n mode in `m_list` and `n_list`
        '''
        if hamada_dphi is None:
            hamada_dphi = self.hamada_dphi

        n_list = np.array(n_list).flatten()
        if len(n_list)!=1 and (len(m_list)!=len(n_list)):
            raise ValueError('Input lists of poloidal and toroidal modes should have the same length, or input only one toroidal mode')
        elif len(n_list) == 1:
            n_list = np.array([n_list[0]]*len(m_list))
        t_array = np.array(range(t_min,t_max+1))
        mode_amplitudes = np.zeros((len(m_list),len(t_array)))
        for t in t_array:
            B = self.get_B_mesh(t)
            if hamada_dphi is None:
                B_n_fft, n_modes, m_modes = self.fft2(B)
            else:
                if abs(hamada_dphi[0] - hamada_dphi[-1]) < 1e-10:
                    hamada_dphi = hamada_dphi[0:-1]
                B_n_fft, n_modes, m_modes = self.fft2(B,hamada_dphi=hamada_dphi)    

            B_n_sorted, n_modes_sorted, m_modes_sorted = self.sort_fft_indices_and_mesh(B_n_fft,n_modes,m_modes)

            i = 0
            for m,n in zip(m_list,n_list):
                m_indices = np.where((m_modes_sorted == m))[0]
                n_indices = np.where((n_modes_sorted == n))[0]
                mode_amplitudes[i][t] = B_n_sorted[m_indices,n_indices].real[0]
                i+=1

        line_list = []
        for j in range(mode_amplitudes.shape[0]):
            line = ax.plot(t_array*dt,mode_amplitudes[j],label=f'{m_list[j]}/{n_list[j]}')
            line_list.append(line)
        ax.set_ylabel('Mode amplitudes (Tesla)')
        ax.set_xlabel('Time (s)')
        ax.set_title('Amplitude of m/n modes in time')
        ax.legend()
        ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
        return line_list

    def field_fourier_amplitude_contour(self,t,m_min,m_max,n_min,n_max,fig,ax,hamada_dphi=None,part='r'):
        '''! Plot the real fourier amplitude of the magnetic field in the fourier space

        @param t The time step during the time evolution
        @param m_min The lower limit of the poloidal harmonics to be included in the contour
        @param m_max The upper limit of the poloidal harmonics to be included in the contour
        @param n_min The lower limit of the toroidal harmonics to be included in the contour
        @param n_max The upper limit of the toroidal harmonics to be included in the contour
        @param fig Matplotlib figure for plotting
        @param ax Matplotlib axis for plotting
        @param hamada_dphi Hamada phase shifts [ntheta]
        @param part The part of the fourier amplitude to be plotted ('r' for real, 'i' for imaginary, 'a' for absolute)
        @result cf The contour plot of the real fourier amplitude of the magnetic field in the fourier space
        @result cbar The colorbar of the contour plot
        '''
        if hamada_dphi is None:
            hamada_dphi = self.hamada_dphi

        if len(np.array(ax).flatten())!=1 or len(np.array(fig).flatten())!=1:
            raise ValueError('For customized plotting, please provide a single figure and a single axis for plotting the a contour.')

        B = self.get_B_mesh(t)
        if hamada_dphi is None:
            B_n_fft, n_modes, m_modes = self.fft2(B)
        else:
            if abs(hamada_dphi[0] - hamada_dphi[-1]) < 1e-10:
                hamada_dphi = hamada_dphi[0:-1]
            B_n_fft, n_modes, m_modes = self.fft2(B,hamada_dphi=hamada_dphi)

        B_n_sorted, n_modes_sorted, m_modes_sorted = self.sort_fft_indices_and_mesh(B_n_fft,n_modes,m_modes)

        n_indices = np.where((n_modes_sorted == n_min) | (n_modes_sorted == n_max))[0]
        m_indices = np.where((m_modes_sorted == m_min) | (m_modes_sorted == m_max))[0]

        n_grid_fft, m_grid_fft = np.meshgrid(n_modes_sorted[n_indices[0]:n_indices[1]+1],m_modes_sorted[m_indices[0]:m_indices[1]+1])

        ax.set_title(f"Real Fourier Amplitude of the Normal Magnetic Field at the Surface at [t] = {t}")
        ax.set_ylabel(r"$n$")
        ax.set_xlabel(r"$m$")
        if part == 'r':
            amplitudes = B_n_sorted[m_indices[0]:m_indices[1]+1,n_indices[0]:n_indices[1]+1].real
            cf = ax.contourf(m_grid_fft,n_grid_fft,amplitudes,levels=50,cmap="viridis")
            cbar = fig.colorbar(cf,label="Mode real amplitude (Tesla)")
        elif part == 'i':
            amplitudes = B_n_sorted[m_indices[0]:m_indices[1]+1,n_indices[0]:n_indices[1]+1].imag
            cf = ax.contourf(m_grid_fft,n_grid_fft,amplitudes,levels=50,cmap="viridis")
            cbar = fig.colorbar(cf,label="Mode imaginary amplitude (Tesla)")
        elif part == 'a':
            amplitudes = abs(B_n_sorted[m_indices[0]:m_indices[1]+1,n_indices[0]:n_indices[1]+1])
            cf = ax.contourf(m_grid_fft,n_grid_fft,amplitudes,levels=50,cmap="viridis")
            cbar = fig.colorbar(cf,label="Mode absolute amplitude (Tesla)")
        cbar.ax.ticklabel_format(style='sci', scilimits=(-3, 3))
        return cf, cbar

    def plot_2D_fourier_amplitude(self,t,harmonics,ax,toroidal_harmonics=True,hamada_dphi=None,x_type='modes',x_mode_min=None,x_mode_max=None,sensor_mesh=None):
        '''! Plot the 2D Fast Fourier Transformed amplitude of the mesh of magnetic values against poloidal/toroidal harmonics/angles

        @param t The time step during the time evolution
        @param harmonics List of (poloidal/toroidal) modes whose amplitudes are to be visualized in y axis [:]
        @param ax Matplotlib axis for plotting
        @param toroidal_harmonics Whether the input `harmonics` is toroidal or poloidal harmonics
        @param hamada_dphi Hamada phase shifts [ntheta]
        @param x_type The variable on x axis ('modes' or 'angles')
        @param x_mode_min The min of the (toroidal/poloidal) mode number to be visualized on x axis
        @param x_mode_max The max of the (toroidal/poloidal) mode number to be visualized on x axis
        @param sensor_mesh Customized sensor signal mesh to be used for FFT and plotting (often from frequency response calculation) [ntheta, nphi]
        @result real_line_list The list of line objects of real amplitudes for each harmonic in `harmonics`
        @result imag_line_list The list of line objects of imaginary amplitudes for each harmonic in `harmonics`
        '''
        if hamada_dphi is None:
            hamada_dphi = self.hamada_dphi

        harmonics = np.array([harmonics]).flatten()
        if x_type not in ['modes','angles']:
            raise ValueError("Unsupported x variable is provided. Accepts 'modes' and 'angles' only.")
        elif x_type == 'modes' and (x_mode_min is None or x_mode_max is None):
            raise ValueError('For x_type == "modes", both mode_min and mode_max should be provided.')
            
        if sensor_mesh is None:
            B = self.get_B_mesh(t)
        elif sensor_mesh.shape != (self.ntheta,self.nphi):
            raise ValueError("The shape of sensor_mesh does not match the sensor dimensions.")
        else:
            B = sensor_mesh
            
        if hamada_dphi is None:
            B_n_fft, n_modes, m_modes = self.fft2(B)
        else:
            if abs(hamada_dphi[0] - hamada_dphi[-1]) < 1e-10:
                hamada_dphi = hamada_dphi[0:-1]
            B_n_fft, n_modes, m_modes = self.fft2(B,hamada_dphi=hamada_dphi)

        import matplotlib.pyplot as plt
        if x_type == 'modes':
            B_n_sorted, n_modes_sorted, m_modes_sorted = self.sort_fft_indices_and_mesh(B_n_fft,n_modes,m_modes)
            cmap = plt.get_cmap('tab10')
            if toroidal_harmonics:
                m_range = np.where((m_modes_sorted == x_mode_min) | (m_modes_sorted == x_mode_max))[0]
                mode_idx = [np.where(n_modes_sorted == harmonic)[0][0] for harmonic in harmonics]
                real_line_list = []
                imag_line_list = []
                for i, idx in enumerate(mode_idx):
                    color = cmap(i/len(mode_idx))
                    rline = ax.plot(m_modes_sorted[m_range[0]:m_range[1]+1],B_n_sorted[m_range[0]:m_range[1]+1,idx].real,color=color,label=f"n={harmonics[i]}, real")
                    iline = ax.plot(m_modes_sorted[m_range[0]:m_range[1]+1],B_n_sorted[m_range[0]:m_range[1]+1,idx].imag,linestyle='--',color=color,label=f"n={harmonics[i]}, imag")
                    real_line_list.append(rline)
                    imag_line_list.append(iline)
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Toroidal Modes at [t] = {t}")
                ax.set_xlabel(r"Poloidal Harmonics ($m$)")
                ax.set_ylabel("Mode Amplitudes (Tesla)")
                ax.set_xticks(range(x_mode_min,x_mode_max+1,2))
                ax.grid()
                ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
                # plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
                return real_line_list, imag_line_list
            else:
                n_range = np.where((n_modes_sorted == x_mode_min) | (n_modes_sorted == x_mode_max))[0]
                mode_idx = [np.where(m_modes_sorted == harmonic)[0][0] for harmonic in harmonics]
                real_line_list = []
                imag_line_list = []
                for i, idx in enumerate(mode_idx):
                    color = cmap(i/len(mode_idx))
                    rline = ax.plot(n_modes_sorted[n_range[0]:n_range[1]+1],B_n_sorted[idx,n_range[0]:n_range[1]+1].real,color=color,label=f"m={harmonics[i]}, real")
                    iline = ax.plot(n_modes_sorted[n_range[0]:n_range[1]+1],B_n_sorted[idx,n_range[0]:n_range[1]+1].imag,linestyle='--',color=color,label=f"m={harmonics[i]}, imag")
                    real_line_list.append(rline)
                    imag_line_list.append(iline)
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Poloidal Modes at [t] = {t}")
                ax.set_xlabel(r"Toroidal Harmonics ($n$)")
                ax.set_ylabel("Mode Amplitudes (Tesla)")
                ax.set_xticks(range(x_mode_min,x_mode_max+1,2))
                ax.grid()
                ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
                # plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
                return real_line_list, imag_line_list
        else:
            if toroidal_harmonics:
                mode_idx = [np.where(n_modes == harmonic)[0][0] for harmonic in harmonics]
                real_line_list = []
                imag_line_list = []
                for i, idx in enumerate(mode_idx):
                    rline = ax.plot(self.theta_list,B_n_fft[:,idx].real,label=f"n={harmonics[i]}, real")
                    iline = ax.plot(self.theta_list,B_n_fft[:,idx].imag,linestyle='--',label=f"n={harmonics[i]}, imag")
                    real_line_list.append(rline)
                    imag_line_list.append(iline)
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Toroidal Modes at [t] = {t}")
                ax.set_xlabel(r"$\theta$ (radians)")
                ax.set_ylabel("Mode Amplitudes (Tesla)")
                ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
                # plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
                return real_line_list, imag_line_list
            else:
                mode_idx = [np.where(m_modes == harmonic)[0][0] for harmonic in harmonics]
                phi_list = np.linspace(0,2*np.pi,self.nphi,endpoint=False)
                real_line_list = []
                imag_line_list = []
                for i, idx in enumerate(mode_idx):
                    if self.helicity == -1:
                        rline = ax.plot(phi_list,np.roll(B_n_fft[idx,:].real[::-1],shift=1),label=f"m={harmonics[i]}, real")
                        iline = ax.plot(phi_list,np.roll(B_n_fft[idx,:].imag[::-1],shift=1),linestyle='--',label=f"m={harmonics[i]}, imag")
                        real_line_list.append(rline)
                        imag_line_list.append(iline)
                    else:
                        rline = ax.plot(phi_list,B_n_fft[idx,:].real,label=f"m={harmonics[i]}, real")
                        iline = ax.plot(phi_list,B_n_fft[idx,:].imag,linestyle='--',label=f"m={harmonics[i]}, imag")
                        real_line_list.append(rline)
                        imag_line_list.append(iline)
                ax.legend()
                ax.set_title(f"2D FFT Amplitudes of Poloidal Modes at [t] = {t}")
                ax.set_xlabel(r"$\phi$ (radians)")
                ax.set_ylabel("Mode Amplitudes (Tesla)")
                ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
                # plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
                return real_line_list, imag_line_list

    def plot_sensor_signal_against_angle(self,t,ax,theta=True):
        '''! Plot the value of normal magnetic fields over theta/phi at at phi/theta = 0

        @param t The time step during the time evolution
        @param ax Matplotlib axis for plotting
        @param theta Plot against theta (True) or phi (False)
        @return line The line object of the plot
        '''

        B_n = self.get_B_mesh(t)
        if theta:
            line = ax.plot(self.theta_list/2/np.pi, B_n[:,0])
            ax.set_title(rf"Magnetic Field on surface @ $\phi$=0 at [t] = {t}")
            ax.set_xlabel(r"$\theta_{norm}$")
            ax.set_ylabel("Magnetic Field (Tesla)")
        else:
            phi_list = np.linspace(0,2*np.pi,self.nphi,endpoint=False)
            line = ax.plot(phi_list/2/np.pi, np.flip(B_n[0,:],axis=1))
            ax.set_title(rf"Magnetic Field on surface @ $\theta$=0 at [t] = {t}")
            ax.set_xlabel(r"$\phi_{norm}$")
            ax.set_ylabel("Magnetic Field (Tesla)")
        ax.ticklabel_format(style='sci', scilimits=(-3,3), axis='y')
        return line

# ---------------------------------------------------------------------------
# Time-dependent drive construction from equilibria
# ---------------------------------------------------------------------------

def make_drive_and_xml_from_eqdsks(time_arr,eqdsk_filenames,xml_filename,drive_filename,eta_values,coil_mask=1,fig=None,ax=None,plot_crosect_time_index=0,plot_ip=False):
    '''! Read in multiple gEQDSK files and create a time-dependent current drive file and an xml file for coil definition.

    @param time_arr Array of time points corresponding to each gEQDSK file
    @param eqdsk_filenames List of paths to gEQDSK files
    @param xml_filename Path to save the xml file defining plasma current filaments
    @param drive_filename Path to save the time-dependent current drive file. Only
           filaments inside the last closed flux surface at one or more time points are
           written (in the same order as the xml `coil_set` entries)
    @param eta_values Resistivities for the conducting mesh, one per mesh region [nregions]
    @param fig Matplotlib figure for plotting
    @param ax Matplotlib axis for plotting
    @param plot_crosect_time_index Index of the time step to plot cross-sections
    @param coil_mask Mask plasma current filament from sensors
    '''
    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon
    from scipy.interpolate import interp1d, RegularGridInterpolator
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from ..TokaMaker.util import read_eqdsk
    if isinstance(time_arr, (int, float)):
        time_arr = [time_arr]
    if isinstance(eqdsk_filenames, str):
        eqdsk_filenames = [eqdsk_filenames]
        
    if len(eqdsk_filenames) != len(time_arr):
        raise ValueError("Number of gEQDSK files does not match number of time points provided.")
    if len(eqdsk_filenames) == 1 and plot_ip:
        raise ValueError("Plotting Ip is only supported when multiple gEQDSK files are provided.")
    if plot_crosect_time_index >= len(time_arr):
        raise ValueError("plot_crosect_time index is out of bounds for the time array.")
    print('Creating xml and driver files...')

    print('Number of drivers: %d' % (len(eqdsk_filenames)))
    I_mesh_list = []
    zero_time = False
    index = 0
    Ip_list = []
    for filename in eqdsk_filenames:
        print("Now reading: " + filename)
        eqdsk_obj = read_eqdsk(filename)
        ip = eqdsk_obj['ip']
        Psi_grid = eqdsk_obj['psirz']

        # 1D arrays from EQDSK
        pprime_1d = eqdsk_obj['pprime']
        ffprim_1d = eqdsk_obj['ffprim'] 
        psi_1d = np.linspace(eqdsk_obj['psimag'],eqdsk_obj['psibry'],eqdsk_obj['nr'])

        # Interpolate to 2D grid
        pprime_interp = interp1d(psi_1d,pprime_1d,bounds_error=False,fill_value="extrapolate")
        ffprim_interp = interp1d(psi_1d,ffprim_1d,bounds_error=False,fill_value="extrapolate")

        PP_grid = pprime_interp(Psi_grid)
        FF_grid = ffprim_interp(Psi_grid)

        # Original R and Z values
        R_vals = np.linspace(eqdsk_obj['rleft'],eqdsk_obj['rleft']+eqdsk_obj['rdim'],eqdsk_obj['nr'])
        Z_vals = np.linspace(eqdsk_obj['zmid']-eqdsk_obj['zdim']/2,eqdsk_obj['zmid']+eqdsk_obj['zdim']/2,eqdsk_obj['nz'])
        # Undersampling to a 60*120 grid
        if not zero_time:
            R_under = np.linspace(eqdsk_obj['rleft'],eqdsk_obj['rleft']+eqdsk_obj['rdim'],60)
            Z_under = np.linspace(eqdsk_obj['zmid']-eqdsk_obj['zdim']/2,eqdsk_obj['zmid']+eqdsk_obj['zdim']/2,60)
            R_grid_under, Z_grid_under = np.meshgrid(R_under, Z_under)
        f_PP = RegularGridInterpolator((Z_vals,R_vals),PP_grid,bounds_error=False,fill_value=np.nan)
        f_FF = RegularGridInterpolator((Z_vals,R_vals),FF_grid,bounds_error=False,fill_value=np.nan)
        PP_grid_under = f_PP((Z_grid_under, R_grid_under))
        FF_grid_under = f_FF((Z_grid_under, R_grid_under))

        r = R_grid_under[0,:]
        z = Z_grid_under[:,0]
        grid_area = (r[1] - r[0]) * (z[1] - z[0])
        I_mesh = (R_grid_under*PP_grid_under+FF_grid_under/(R_grid_under*mu0))*grid_area

        # Determine which points are within the last closed flux surface
        polygon = Polygon(eqdsk_obj['rzout'])

        # Trim because we only care about inside the last closed flux 
        insideFlux = np.zeros((z.shape[0], r.shape[0]))
        coord_in = [[],[]]
        for i in range(R_grid_under.shape[1]):
            for j in range(Z_grid_under.shape[0]):
                R_tmp = r[i]
                Z_tmp = z[j]
                point = Point(R_tmp,Z_tmp)
                insideFlux[j,i] = polygon.contains(point)
                if(insideFlux[j,i] == 1):
                    coord_in[0].append([R_tmp])
                    coord_in[1].append([Z_tmp])
                else:
                    I_mesh[j,i] = 0.0
        coord_in = np.array(coord_in)
        coord_in = np.squeeze(coord_in)
        print('The length of the trimmed data is: %d' % (coord_in[0].shape[0]))

    
        Ip_tot = np.sum(I_mesh)
        Ip_list.append(Ip_tot)
        if not zero_time:
            Ip_list.append(Ip_tot)
        print(f'Ip from eqdsk: {ip}; Ip from sum of currents on grids: {Ip_tot}')

        I_mesh_list.append(I_mesh.copy())
        zero_time = True
        if index == plot_crosect_time_index:
            if ax is not None and fig is not None:
                def _edges_from_centers(arr1d):
                    """Convert 1D cell centers to cell-edge coordinates (len = N+1)."""
                    arr1d = np.asarray(arr1d)
                    edges = np.empty(arr1d.size+1,dtype=float)
                    # interior edges
                    edges[1:-1] = 0.5*(arr1d[:-1]+arr1d[1:])
                    # extrapolated outer edges
                    edges[0] = arr1d[0]-0.5*(arr1d[1]-arr1d[0])
                    edges[-1] = arr1d[-1]+0.5*(arr1d[-1]-arr1d[-2])
                    return edges
                R_edges = _edges_from_centers(R_under)
                Z_edges = _edges_from_centers(Z_under)
                I_mesh_kA = I_mesh * 1e-3
                vmin = np.nanmin(I_mesh_kA)
                vmax = np.nanmax(I_mesh_kA)
                norm = Normalize(vmin=vmin, vmax=vmax)
                pcm = ax.pcolormesh(R_edges,Z_edges,I_mesh_kA,cmap="viridis",norm=norm,shading="flat",edgecolors='black',linewidth=0.6)
                ax.set_xlabel("R (m)")
                ax.set_ylabel("Z (m)")
                ax.set_aspect("equal",adjustable="box")
                ax.grid(True)
                ax.set_axisbelow(True)
                cbar = fig.colorbar(pcm,ax=ax,pad=0.02,label="Flat-top filament current (kA)")
            elif (ax is not None and fig is None) or (ax is None and fig is not None):
                raise ValueError("Both ax and fig must be provided for plotting, or neither.")
        index += 1

    # Only keep filaments that carry current inside the LCFS at any time point, so the
    # model does not include (60*60 - kept) zero-current coil_sets
    active = np.zeros(R_grid_under.shape, dtype=bool)
    for I_mesh in I_mesh_list:
        active |= (I_mesh != 0.0)
    n_active = int(np.count_nonzero(active))
    print('Retaining %d of %d filaments inside the LCFS across all time points' % (n_active, active.size))

    # Write the xml file (coil_set order defines the drive column order)
    fid = open(xml_filename,'w')
    fid.write('<oft>\n')
    fid.write('  <thincurr>\n')
    fid.write('    <eta>%s</eta>\n' % ', '.join('%.3e' % eta for eta in eta_values))
    fid.write('    <icoils>\n')
    for i in range(Z_grid_under.shape[0]):
        for j in range(R_grid_under.shape[1]):
            if not active[i,j]:
                continue
            fid.write('      <coil_set sens_mask="%d">\n' % (coil_mask))
            fid.write('        <coil>%2.16f, %2.16f</coil>\n' % (R_grid_under[i,j], Z_grid_under[i,j]))
            fid.write('      </coil_set>\n')
    fid.write('    </icoils>\n')
    fid.write('  </thincurr>\n')
    fid.write('</oft>\n')
    fid.close()

    # Write the drive file: a hold row at t=0 followed by one row per time point
    fid0 = open(drive_filename, 'w')
    # ThinCurr expects 'ncols ntimes' (thincurr_td.F90: READ ncols,ntimes then ALLOCATE(ntimes,ncols))
    fid0.write('%d %d\n' % (n_active+1,len(eqdsk_filenames)+1))
    for row_time, I_mesh in zip([0.0]+list(time_arr), [I_mesh_list[0]]+I_mesh_list):
        fid0.write('%f' % (row_time))
        for i in range(Z_grid_under.shape[0]):
            for j in range(R_grid_under.shape[1]):
                if active[i,j]:
                    fid0.write(' %f' % (I_mesh[i,j]))
        fid0.write('\n')
    fid0.close()
    if plot_ip:
        plt.figure(figsize=(10, 8))
        tarr = [0.0]+time_arr
        plt.plot(tarr,np.array(Ip_list)/1e6)
        plt.xlabel('Time (s)')
        plt.ylabel('Plasma Current (MA)')
    if ax is not None and fig is not None:
        return pcm, cbar
    

def drive_to_array(drive_filename):
    '''! Read the time-dependent current drive file and return it as a numpy array.
    
    @param drive_filename Path to the time-dependent current drive file
    @result A numpy array where each row corresponds to a time step array; the first element is the time and the subsequent elements are the currents for each coil_set.
    '''
    import re
    
    with open(drive_filename, 'r') as f:
        first_line = f.readline().strip()
        match = re.match(r'(\d+)\s+(\d+)', first_line)
        if match:
            ncols, nrows = map(int, match.groups())
        else:
            raise ValueError("Invalid format in drive file header")
        return np.array([list(map(float, line.split())) for line in f])


def triangular_waveform(twidth,amplitude,curr_arr,ax=None):
    '''! Build a three-point triangular current waveform from a steady-state current vector

    The waveform ramps every coil current linearly from its steady value up to
    `(1+amplitude)` times that value at `twidth/2`, then back down at `twidth`.

    @param twidth Full width of the triangular pulse in seconds
    @param amplitude Fractional current increase at the peak (|amplitude| <= 1)
    @param curr_arr A single row of a drive waveform: time in element 0, currents after
    @param ax Matplotlib axis to plot the waveform on, or `None` to skip plotting
    @result `tri_waveform` Waveform array [3,ncols] suitable for `ThinCurr.run_td(coil_currs=...)`
    '''
    if np.abs(amplitude) > 1:
        raise ValueError("Amplitude must be less than or equal to 1")
    curr_arr_mid = curr_arr * (amplitude+1)
    curr_arr_end = np.copy(curr_arr)
    curr_arr_start = np.copy(curr_arr)

    curr_arr_end[0] = twidth
    curr_arr_mid[0] = twidth/2
    curr_arr_start[0] = 0.0
    tri_waveform = np.array([curr_arr_start, curr_arr_mid, curr_arr_end])
    yvals = np.array([0,amplitude,0])
    if ax is not None:
        from matplotlib.ticker import PercentFormatter
        ax.plot(tri_waveform[:,0]*1e3, yvals, marker='o')
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel(r'% increase of $I_{p}$')
        ax.set_title('Triangular Waveform of Current')
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=1))
        data_min, data_max = float(yvals.min()), float(yvals.max())
        pad = max(1e-3, 0.02 * (data_max - data_min))
        ax.set_ylim(data_min - pad, data_max + pad)
        ticks = list(ax.get_yticks())
        if not any(abs(t_ - data_max) < 1e-12 for t_ in ticks):
            ticks.append(data_max)
            ticks = np.unique(np.array(ticks))
            ticks.sort()
            ax.set_yticks(ticks)
    return tri_waveform

# ---------------------------------------------------------------------------
# ThinCurr model preparation helpers
# ---------------------------------------------------------------------------

def shift_coils_in_xml(filename, output_filename, dx=0.0, dy=0.0, dz=0.0):
    '''! Shift the coordinates of every `<coil>` element in a ThinCurr XML file

    Strict: raises `ValueError` if any `<coil>` lacks an `npts` attribute or has a
    line that is not exactly three comma-separated coordinates. Nothing is written
    unless every coil passes.

    @param filename ThinCurr XML file to read
    @param output_filename Path to write the shifted XML to
    @param dx Shift applied to the X coordinate of every point [m]
    @param dy Shift applied to the Y coordinate of every point [m]
    @param dz Shift applied to the Z coordinate of every point [m]
    @result `(coils_shifted, points_shifted)` Number of coils and coordinate points shifted
    '''
    print(f"Reading: {filename}")
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    total_points_shifted = 0
    coils_shifted = 0
    
    # Find all <coil> tags anywhere in the tree
    for coil_idx, coil in enumerate(root.iter('coil')):
        if coil.text is None or not coil.text.strip():
            raise ValueError(f"Error in coil #{coil_idx + 1}: Coil block is completely empty.")
        
        # STRICT CHECK 1: Must have 'npts' attribute
        if 'npts' not in coil.attrib:
            raise ValueError(f"Error in coil #{coil_idx + 1}: Missing 'npts' attribute. Cannot verify if 3D.")
            
        raw_lines = coil.text.strip().split('\n')
        shifted_lines = []
        
        for line_idx, line in enumerate(raw_lines):
            if not line.strip():
                continue
                
            parts = line.split(',')
            
            # STRICT CHECK 2: Must have exactly 3 coordinates (X, Y, Z)
            if len(parts) != 3:
                raise ValueError(
                    f"Error in coil #{coil_idx + 1}, line {line_idx + 1}: "
                    f"Expected 3 coordinates (X, Y, Z), but found {len(parts)}.\n"
                    f"Line content: '{line.strip()}'"
                )
            
            # Apply the shifts
            try:
                x = float(parts[0].strip()) + dx
                y = float(parts[1].strip()) + dy
                z = float(parts[2].strip()) + dz
            except ValueError:
                raise ValueError(f"Error in coil #{coil_idx + 1}: Non-numeric coordinate found on line {line_idx + 1}.")
            
            shifted_lines.append(f"{x}, {y}, {z}")
            total_points_shifted += 1
        
        # Update the coil text
        coil.text = "\n" + "\n".join(shifted_lines) + "\n"
        coils_shifted += 1
        
    # If we made it through EVERY coil without crashing, write the file
    tree.write(output_filename, encoding='utf-8', xml_declaration=False)
    
    print("\n--- Strict Shift SUCCESS ---")
    print(f"Applied shifts -> dx: {dx}m, dy: {dy}m, dz: {dz}m")
    print(f"Total Coils shifted: {coils_shifted}")
    print(f"Total individual coordinates shifted: {total_points_shifted}")
    print(f"Saved to: {output_filename}\n")
    return coils_shifted, total_points_shifted
    


def parse_coils_xml(filename):
    '''! Read the coil coordinates from a ThinCurr XML file

    @param filename ThinCurr XML file to read
    @result `coils` List of numpy coordinate arrays, one per `<coil>` element [npts,3]
    '''
    tree = ET.parse(filename)
    root = tree.getroot()
    coils = []
    
    # Iterate through all <coil> elements in the XML
    for coil in root.findall('.//coil'):
        # The coordinates are raw text, separated by newlines
        text = coil.text.strip()
        points = []
        for line in text.split('\n'):
            line = line.strip()
            if line: # Ensure the line isn't empty
                coords = [float(val) for val in line.split(',')]
                points.append(coords)
        coils.append(np.array(points))
    return coils




def find_coil_current_column(coil_name,header_to_idx,header_template='{name} current [a/turn]',
                             ignore_suffixes=('_feed',)):
    '''! Find the CSV column holding the current waveform for a named coil set

    Coil-set names in a ThinCurr model often differ from the CSV headers by feed
    markers or per-loop numbering, so candidates are tried in order: the name as
    given, the name with each of `ignore_suffixes` removed, and each of those with
    a trailing `_<number>` removed. All matching is case-insensitive.

    @param coil_name Coil-set name to look up
    @param header_to_idx Mapping of lower-cased CSV header to column index
    @param header_template Header pattern, with `{name}` replaced by each candidate name
    @param ignore_suffixes Name fragments to strip when the exact name has no match
    @result Column index of the first matching header, or `None` if no candidate matches
    '''
    name_lower = coil_name.lower()
    candidates = [name_lower]
    stripped = name_lower
    for suffix in ignore_suffixes:
        stripped = stripped.replace(suffix.lower(), '')
    if stripped != name_lower:
        candidates.append(stripped)
    for base in list(candidates):
        no_number = re.sub(r'_\d+$', '', base)
        if no_number not in candidates:
            candidates.append(no_number)
    for candidate in candidates:
        target = header_template.format(name=candidate)
        if target in header_to_idx:
            return header_to_idx[target]
    return None


def add_gpec_coils_to_xml(xml_filename,nc_filename,output_filename,prefixes=None,
                          keep_existing=True,max_npts=None,sens_mask='0',verbose=True):
    '''! Append coil geometry from a GPEC coil-output file to a ThinCurr XML model

    Each GPEC coil group appears as `<prefix>_x`, `<prefix>_y`, `<prefix>_z` variables of
    shape [npts,nsubsystems,ncoils]. Groups whose `<prefix>_current` holds more than one
    value are split into one `coil_set` per loop, so that each can be driven independently;
    all others share a single `coil_set`. Loops that are entirely NaN or zero (padding in
    the GPEC file) are skipped, and a group that contributes no valid loop is removed again.

    @param xml_filename ThinCurr XML model to read
    @param nc_filename GPEC coil-output netCDF file to read geometry from
    @param output_filename Path to write the augmented ThinCurr XML to
    @param prefixes Iterable of GPEC coil-group prefixes to include, or `None` for all
    @param keep_existing Keep the `icoils` already in the model (`False` clears them first,
           e.g. to write a coils-only model without the plasma filaments)
    @param max_npts If given, loops longer than this are downsampled, always keeping the first and last point
    @param sens_mask Value of the `sens_mask` attribute written on each `coil_set`
    @param verbose Print a summary of the coil sets that were added
    @result `added_coil_sets` Names of the coil sets added, in the order ThinCurr will
            assign waveform columns to them
    '''
    from scipy.io import netcdf_file

    tree = ET.parse(xml_filename)
    root = tree.getroot()
    thincurr = root.find('./thincurr')
    if thincurr is None:
        raise ValueError("Could not find <thincurr> in %s" % xml_filename)

    icoils = thincurr.find('./icoils')
    if icoils is None:
        icoils = ET.SubElement(thincurr, 'icoils')
    elif not keep_existing:
        for child in list(icoils):
            icoils.remove(child)

    ds = netcdf_file(nc_filename, 'r', mmap=False)
    try:
        wanted = None if prefixes is None else set(prefixes)
        coil_prefixes = set()
        for var_name in ds.variables:
            var_str = var_name.decode('utf-8') if isinstance(var_name, bytes) else var_name
            if not var_str.endswith('_x'):
                continue
            prefix = var_str[:-2]
            if wanted is None or prefix in wanted:
                coil_prefixes.add(prefix)
        if wanted is not None:
            missing = wanted - coil_prefixes
            if missing:
                raise ValueError("Coil prefixes not present in %s: %s"
                                 % (nc_filename, sorted(missing)))
        if verbose:
            print("Found %d coil groups in %s" % (len(coil_prefixes), nc_filename))

        valid_coils_count = 0
        added_coil_sets = []
        for prefix in sorted(coil_prefixes):
            x_data = ds.variables[prefix+'_x'][...]
            y_data = ds.variables[prefix+'_y'][...]
            z_data = ds.variables[prefix+'_z'][...]
            nw_value = float(ds.variables[prefix+'_nw'][...])*-1.0
            n_points, n_subsystems, n_coils = x_data.shape

            current_var = prefix+'_current'
            split_coils = False
            if current_var in ds.variables:
                split_coils = ds.variables[current_var][...].size > 1

            coil_set = None
            if not split_coils:
                coil_set = ET.SubElement(icoils, 'coil_set',
                                         {'sens_mask': sens_mask, 'name': prefix})
            group_valid_loops = 0
            split_coil_counter = 1
            for sub in range(n_subsystems):
                for c in range(n_coils):
                    coil_x = x_data[:, sub, c]
                    coil_y = y_data[:, sub, c]
                    coil_z = z_data[:, sub, c]
                    if max_npts is not None and n_points > max_npts:
                        # always keep the first and last point
                        indices = np.linspace(0, n_points-1, max_npts, dtype=int)
                        coil_x = coil_x[indices]
                        coil_y = coil_y[indices]
                        coil_z = coil_z[indices]
                        actual_npts = max_npts
                    else:
                        actual_npts = n_points
                    if np.all(np.isnan(coil_x)) or np.all(coil_x == 0):
                        continue

                    group_valid_loops += 1
                    valid_coils_count += 1
                    if split_coils:
                        individual_name = "%s_%d" % (prefix, split_coil_counter)
                        current_coil_set = ET.SubElement(
                            icoils, 'coil_set', {'sens_mask': sens_mask, 'name': individual_name})
                        added_coil_sets.append(individual_name)
                        split_coil_counter += 1
                    else:
                        current_coil_set = coil_set

                    coil_element = ET.SubElement(current_coil_set, 'coil')
                    coil_element.set('npts', str(actual_npts))
                    coil_element.set('scale', str(nw_value))
                    coord_lines = [f"{x}, {y}, {z}"
                                   for x, y, z in zip(coil_x, coil_y, coil_z)]
                    coil_element.text = "\n"+"\n".join(coord_lines)+"\n"

            if not split_coils:
                if group_valid_loops == 0:
                    icoils.remove(coil_set)
                else:
                    added_coil_sets.append(prefix)
    finally:
        ds.close()

    if hasattr(ET, 'indent'):
        ET.indent(tree, space="  ", level=0)
    tree.write(output_filename, encoding='utf-8', xml_declaration=False)

    if verbose:
        print("Processed %d coils into %s" % (valid_coils_count, output_filename))
        for i, name in enumerate(added_coil_sets, start=1):
            print("%03d. %s" % (i, name))
    return added_coil_sets


def append_coil_currents_to_drive(drive_filename,csv_filename,output_filename,coil_names,
                                  fill_value='0.000000',header_template='{name} current [a/turn]',
                                  ignore_suffixes=('_feed',),verbose=True):
    '''! Append machine coil-current waveforms from a CSV to an existing `.drive` file

    Column order follows `coil_names`, which must be the order the coil sets appear in the
    ThinCurr XML (as returned by `add_gpec_coils_to_xml`), because that is the order
    ThinCurr assigns waveform columns. Coil names are matched to CSV headers by
    `find_coil_current_column`; a name with no match, or a time row absent from the CSV, is filled
    with `fill_value` so the matrix stays rectangular.

    @param drive_filename Existing `.drive` file to extend
    @param csv_filename CSV of coil currents, one column per coil named `<coil> current [A/turn]`
    @param output_filename Path to write the extended `.drive` file to
    @param coil_names Coil-set names to append, in ThinCurr column order
    @param fill_value Value written where a coil or time row has no CSV entry
    @param header_template Header pattern used to match coils to CSV columns (see `find_coil_current_column`)
    @param ignore_suffixes Name fragments ignored during matching (see `find_coil_current_column`)
    @param verbose Print a summary and warn about unmatched coils
    @result `missing_columns` Names that had no CSV match and were filled with `fill_value`
    '''
    import csv as _csv

    with open(drive_filename, 'r') as fid:
        drive_lines = fid.readlines()
    header_dims = drive_lines[0].strip().split()
    data_line_1 = drive_lines[1].strip().split()
    if len(data_line_1) == int(header_dims[0]):
        old_num_cols = int(header_dims[0])
        num_rows = int(header_dims[1])
        cols_first = True
    else:
        num_rows = int(header_dims[0])
        old_num_cols = int(header_dims[1])
        cols_first = False
    drive_matrix = [line.strip().split() for line in drive_lines[1:] if line.strip()]

    with open(csv_filename, 'r', encoding='utf-8-sig') as fid:
        reader = _csv.reader(fid)
        csv_headers = next(reader)
        csv_data = list(reader)
    header_to_idx = {h.strip().lower(): i for i, h in enumerate(csv_headers)}

    new_columns_data = [[] for _ in range(num_rows)]
    missing_columns = set()
    for coil_name in coil_names:
        col_idx = find_coil_current_column(coil_name, header_to_idx)
        if col_idx is None:
            missing_columns.add(coil_name)
        for row_idx in range(num_rows):
            if col_idx is not None and row_idx < len(csv_data):
                new_columns_data[row_idx].append(csv_data[row_idx][col_idx])
            else:
                new_columns_data[row_idx].append(fill_value)

    new_total_cols = old_num_cols+len(coil_names)
    with open(output_filename, 'w') as fid:
        if cols_first:
            fid.write("%d %d\n" % (new_total_cols, num_rows))
        else:
            fid.write("%d %d\n" % (num_rows, new_total_cols))
        for r in range(num_rows):
            fid.write("%s %s\n" % (" ".join(drive_matrix[r]), " ".join(new_columns_data[r])))

    if verbose:
        if missing_columns:
            print("WARNING: no CSV column for: %s" % ", ".join(sorted(missing_columns)))
            print("These were filled with %s to keep the matrix rectangular." % fill_value)
        print("Appended %d columns -> %s (%d cols, %d rows)"
              % (len(coil_names), output_filename, new_total_cols, num_rows))
    return missing_columns
