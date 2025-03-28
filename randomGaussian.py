import os
import logging
from logging import WARNING, ERROR, CRITICAL, INFO
from typing import Optional
import numpy as np
from sys import stdout

logger = logging.getLogger("RandomDistributions")
handler   = logging.StreamHandler(stdout)
logger.setLevel(INFO)#DEBUG)
formatter = logging.Formatter('%(asctime)s - [%(levelname)s]: %(message)s')

handler   = logging.StreamHandler(stdout)
handler.setFormatter(formatter)
logger.addHandler(handler)

np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
matplotlib    = True
scipy         = True
cythonFortran = True
papersettings = True

try: 
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from matplotlib.ticker import ScalarFormatter
    from matplotlib.gridspec import GridSpec 
    try:
        import paper_settings as pg 
        paper = pg.PaperSettings()
        plt.rcParams.update({
                "text.usetex": True, **paper.mplfigfontstyle
                #"font.family": "Helvetica"
                })
    except:
        logger.info("No paper settings found.") 
        papersettings = False
    
except: 
    matplotlib = False
    logger.info("No matplotlib package found.")
try:
    from scipy.fftpack import fftfreq, fftn, ifftn
except:
    fftfreq, fftn, ifftn = np.fft.fftfreq, np.fft.fftn, np.fft.ifftn
    scipy = False
    logger.info("No scipyfft found, using numpy.fft.")
try: 
    from cython_fortran_file import FortranFile

except:
    from scipy.io import FortranFile
    FortranFile.write_vector = FortranFile.write_record
    cythonFortran = False
    logger.info("No cython_fortran_file found, using scipy.io.")




themes = {"div" : dict(cmap = "RdBu", color="black"),
          "dark": dict(cmap = "cividis", color="white")}
theme    = themes["div"]
txt_color = theme["color"]
cmap      = theme["cmap"]



class RandomDistribution:
    logger = logging.getLogger("RandomDistributions")
    def __init__(self, 
                 ndim: int,
                 size: int|tuple|np.ndarray,
                 mean: float|int  = 0,
                 rms : float|int  = 0.1,
                 Boxlen: float|int= 1,
                 seed  : int = 0,
                 kmin :Optional[float|int] = None,
                 kmax :Optional[float|int] = None,
                 lmin :Optional[float|int] = None,
                 lmax :Optional[float|int] = None,
                 nkbins = 100
                 ):
        self.ndim    = ndim 
        self.size    = self._toarray(size, ndim)
        assert ndim == len(self.size)
        self.mean    = mean
        self.rms     = rms
        self.Boxlen  = self._toarray(Boxlen, ndim)
        assert ndim == len(self.Boxlen)
        
        self.dx      = Boxlen/self.size
        self.seed    = seed
        
        self.kmin    = kmin
        self.kmax    = kmax
        self.lmin    = lmin
        self.lmax    = lmax
         
        if self.lmax is not None: 
            if self.kmin is not None: logger.warning("supplied kmin will be ignored if you give lmax.")
            self.kmin = 2.*np.pi/self.lmax
        if self.lmin is not None: 
            if self.kmax is not None: logger.warning("supplied kmax will be ignored if you give lmin.")
            self.kmax = 2.*np.pi/self.lmin
        pass
    
    @staticmethod
    def _toarray(arg, ndim:int):
        if not isinstance(arg, np.ndarray):
            if not isinstance(arg, list): 
                return np.array([arg]*ndim) 
            else: 
                
                return np.array(arg)
        else:
            return arg
        
    def distribution(self, pk_function = lambda k : k**(-5./3.), nbins = 100, kcorr = None,sigmak = 1):
           
        # create the frequencies from the signal's number with np.fttfreq 

        # this is equivalent to     ki = [kx = np.fft.fftfreq(nx),...,   kn = np.fft.fftfreq(nn)]
        ndim = self.ndim 
        size = self.size
        seed = self.seed
        kmin, kmax = self.kmin, self.kmax
        ki = [  fftfreq(n, d = self.dx[i]).reshape([-1 if i == j else 1 for j in range(ndim)]) for i, n in enumerate(size) ]
        # create agrid dynamically, similar to
        # kkx,..., kkn = np.meshgrid(kx,..., kn, indexing="ij")
        #k2 = kkx**2 +...+ kkn**2
        # create the k grid
        k2 = sum(k**2 for k in ki)
        k  = np.sqrt(k2)
        # we set the first element to 1 to avoid divergent values in case of powerlaw specturms
        # this is equivalent to k[0_1, ..., 0_n] = 1.0
        k.flat[0] = 1.0

        # now we generate our power-specturm 
        np.random.seed(seed)
        if kmin is None and kmax is None:
            pk = pk_function(k)#k**alpha
            if kcorr is not None:
                Window = np.exp(-((k - kcorr) / sigmak )**2 )
                pk *= Window
                
            # now we normalize the specturm base on the grid volume     
            pk /= np.sum(pk)     # over the sum 
            pk *= np.prod(size)  # and scale by total number of grid points

            # Now that we have the specturm, we now want to randomize the complex amplutudes of the signal. 
            # therefore, we generate a gaussian noise
                
            noise = np.random.normal(0, 1, size) + 1j * np.random.normal(0, 1, size)
            # we now randomize our nodes and thake the amplitudes, nowing that
            # |A(k)|**2 = P(k) -> |A(k)| = (P(k))**1/2 

        else:
        ### do computations on a specific log-spaced range of k
        # mmmh....
        
        
            kbins = np.logspace(np.log10(kmin), np.log10(kmax),nbins)
            pkbins = pk_function(kbins)

            # Assign bin indices to k values
            k_indices = np.digitize(k.flat, kbins, right=True)

            # Initialize power spectrum
            pk = np.zeros_like(k)

            # Assign power spectrum to bins
            for i, pkbin in enumerate(pkbins):
                pk.flat[k_indices == i] = pkbin
            if kcorr is not None:
                Window = np.exp(-((k - kcorr) / sigmak )**2 )
                pk *= Window
        
            # Normalize the power spectrum
            pk /= np.sum(pk)
            pk *= np.prod(size)

            # Generate complex Gaussian noise
            noise = np.random.normal(0, 1, size) + 1j * np.random.normal(0, 1, size)
            #noise[k>kmax]=0.0
            # Apply the power spectrum
        
        #filter_function = self.top_hat_filter(k, kmin, kmax)
        #pk*=filter_function
        fourier_field = noise * np.sqrt(pk)
        
        # back to real space with the inverse of the fourier distribution
        # but we only keep the real component 
        
        field = ifftn(fourier_field).real  
        
        # normalize again the field 
        field = field / np.std(field)
        
        # rescale the field 
        field *= self.rms
        
        # shift to the desired mean 
        field += self.mean
        return field
    @staticmethod
    def top_hat_filter( k, k_min, k_max, delta_k = None):
        """
        Smooth top-hat filter for the power spectrum.
        
        Parameters:
            k (numpy.ndarray): Array of wavenumbers.
            k_min (float): Minimum wavenumber for the filter.
            k_max (float): Maximum wavenumber for the filter.
            delta_k (float): Smoothing width of the filter transition.

        Returns:
            numpy.ndarray: Filter values at each wavenumber.
        """
        # Smooth transitions using hyperbolic tangent
        if delta_k is None: delta_k = 0.8*(k_max-k_min)
        low_pass = 0.5 * (1 + np.tanh((k - k_min) / delta_k))
        high_pass = 0.5 * (1 - np.tanh((k - k_max) / delta_k))
        return low_pass * high_pass
    
    @staticmethod
    def compute_pk(field, kmax = None, kmin = None):
        
        size = field.shape
        gridvolume = np.prod(size)
    
        # compute the fourier transform
        fourier_field = np.abs(fftn(field))**2  # Take the squared magnitude

        ki = [  fftfreq(n).reshape([-1 if i == j else 1 for j in range(ndim)]) for i, n in enumerate(size) ]
        # create agrid dynamically, similar to
        # kkx,..., kkn = np.meshgrid(kx,..., kn, indexing="ij")
        #k2 = kkx**2 +...+ kkn**2
        # create a k grid
        k2 = sum(k**2 for k in ki)
        k  = np.sqrt(k2)
        #k.flat[0] = 0.0
        
        
        # create radial kbins log spaced 
        kflat     = k.flatten()
        pkflat    = fourier_field.flatten()
        if kmin is None: kmin  = min((1./n for n in size))/2 #np.min(kflat[kflat > 0])
        if kmax is None: kmax  = k.max()*2
    
        kbins = 10**np.linspace(np.log10(kmin), np.log10(kmax), 20 )
        
        # Assign each kflat to a bin
        power_sum, kbins     = np.histogram(kflat, bins=kbins, weights=pkflat)
        count_per_bin, _     = np.histogram(kflat, bins=kbins)

        # Avoid division by zero by using np.divide with the `where` parameter
        power_spectrum = np.divide(power_sum, count_per_bin, out = np.zeros_like(power_sum), where = count_per_bin > 0)
        # normalize power specturm 
        power_spectrum /= gridvolume
        pkflat          /= gridvolume
        # Compute the bin centers for plotting
        bin_centers = 0.5 * (kbins[:-1] + kbins[1:])
        valid = count_per_bin > 0
        bin_centers    = bin_centers[valid]
        power_spectrum = power_spectrum[valid]
        
        valid = kflat < kmax
        kflat    = kflat[valid]
        pkflat   = pkflat[valid]
        valid = kflat > kmin
        kflat    = kflat[valid]
        pkflat   = pkflat[valid]
        
        return bin_centers, power_spectrum, kflat, pkflat
    def save(self, field, path="./data", name = "RandomField.dat"):
        if not os.path.exists(path): os.makedirs(path)
        path += "/"
        file = FortranFile(path + name, mode="w")
        file.write_vector((np.array([self.ndim]).flatten()), "i")
        file.write_vector(np.array([self.size], "i").flatten())
        #for line in np.flip(gaussian_field.T,axis=1):
        
        if self.ndim==2:
            for line in np.rot90(field.T):
                file.write_vector(line)
        elif self.ndim==3:
            mygauss = field.copy(order = "F")# np.rot90(gaussian_field, axes=(0,2))
            file.write_vector(mygauss.T)
            return 
        
    def save_vtk(self, field, path= "./data", name = "RandomFieldVTK"):
        from pyevtk.hl import gridToVTK
        
        if not os.path.exists(path): os.makedirs(path)
        path += "/"
        
        ndim  = self.ndim
        size  = self.size
        dx    = self.dx
       
        coord = [np.arange(size[i])*dx[i] - 0.5 * dx[i] for i in range(ndim)]
        
        gridToVTK(
                path + name,
                coord[0], coord[1], coord[2],
                pointData={"field": gaussian_field}
            )
        return   
        
    def plot(self, field, path = "./plots", triplot=True, 
            pk_function = lambda k, k0: (k/k0)**(-5./3.), fig = None,single_cbar=False,
            **kwargs
            ):
        if not os.path.exists(path): os.makedirs(path)
        path +="/"
        first = True
        ndim  = self.ndim
        size  = self.size
        dx    = self.dx
        gaussian_field
        kmag, pk, kflat, pkflat = self.compute_pk(field)
        
        coord = [np.arange(size[i])*dx[i] - 0.5 * dx[i] for i in range(ndim)]
        x_coord = coord[1]  # Second axis (columns)
        # Y-coordinates should correspond to the first axis of `field`
        y_coord = coord[0]
        
        if first:
            if fig is  None:
                fig      = plt.figure(plt.gcf().number+1)
            gs       = GridSpec(3, 3)
            ax1      = fig.add_subplot(gs[1:3, :2])
            ax_xDist = fig.add_subplot(gs[0,   :2] )
            cbaxes = inset_axes( ax1, width = "70%" , height      = "5%", loc = "upper center")    
            
            if papersettings: fig.set_size_inches(2*paper.onecol,2*paper.onecol)
            
            ax_yDist = fig.add_subplot(gs[1:3, 2])
            first = False
        else: 
            for ax in fig.axes: ax.cla()
        if ndim > 1:
             
            if ndim == 2:              
                h = ax1.pcolormesh(x_coord, y_coord, field, cmap = cmap,shading="gouraud", **kwargs)#, vmin=vmin, vmax=vmax)
            else:        
                cut = size[2]//2                
                toplot = field[:, :, cut]
                ext = [toplot.min(),toplot.max()]
                vmin,vmax = ext
                x_coord, y_coord = coord[1], coord[0]
        
                h = ax1.pcolormesh(x_coord, y_coord,toplot , cmap = cmap, vmin=vmin, vmax=vmax,shading="gouraud",**kwargs)
                if triplot:
                    if single_cbar: ext = [(field.flatten()).min(),(field.flatten()).max()]
                    
                    fig22,axx = plt.subplots(2,3, figsize=(12, 4), num = plt.gcf().number + 2)
                    #[ [ax11, ax12], [ax21 , ax22]]  = axx
                    cut = size[2]//2  
                    [ ax11, ax12, ax21,axh11, axh12, axh21  ]  = axx.flatten()
                    toplot = field[:, :, cut]
                    if not single_cbar: ext = [toplot.min(),toplot.max()]
                    x_coord, y_coord = coord[1], coord[0]
                    ax11.pcolormesh(x_coord, y_coord, toplot , cmap = cmap, vmin=vmin, vmax=vmax,shading="gouraud",**kwargs)
                    ax11.set_xlabel("x")
                    ax11.set_ylabel("y")
                    ax11.set_aspect("equal")  
                    axh11.hist(toplot.flat, log=True, bins=50, color = "red" )
                    #print(toplot.flatten().mean(),toplot.flatten().std())
                    #toplot = np.sum(gaussian_field, axis=1)#[:, :, cut]
                    cut = size[1]//2  
                    toplot = field[:, cut, :]
                    x_coord, y_coord = coord[2], coord[0]
                    #toplot = gaussian_field[:, :, cut+1]
                    if not single_cbar: ext = [toplot.flatten().min(),toplot.flatten().max()]
                    vmin,vmax = ext
                    
                    ax12.pcolormesh(x_coord, y_coord,toplot , cmap = cmap, vmin=vmin, vmax=vmax,shading="gouraud",**kwargs)
                    ax12.set_xlabel("x")
                    ax12.set_ylabel("z")
                    ax12.set_aspect("equal")  
                    axh12.hist(toplot.flat, log=True, bins=50, color = "red" )
                    #print(toplot.flatten().mean(),toplot.flatten().std())
                    #toplot = np.sum(gaussian_field, axis=2)#[:, :, cut]
                    cut = size[0]//2  
                    toplot = field[cut, :, :]
                    x_coord, y_coord = coord[2], coord[1]
                    if not single_cbar: ext = [toplot.min(),toplot.max()]
                    vmin,vmax = ext                        
                    ax21.pcolormesh(x_coord, y_coord,toplot , cmap = cmap, vmin=vmin, vmax=vmax,shading="gouraud",**kwargs)
                    ax21.set_xlabel("y")
                    ax21.set_ylabel("z")
                    ax21.set_aspect("equal")  
                    axh21.hist(toplot.flat, log=True, bins=50, color = "red" )
                    #print(toplot.flatten().mean(),toplot.flatten().std())
                    fig22.tight_layout()
                    fig22.savefig(path+"triplot_%05d"%( self.seed))
            cb 	   = fig.colorbar( h, cax   = cbaxes, orientation = "horizontal")
            cb.set_label(label = r"$\delta(x,y)$", color = txt_color)
            plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color = txt_color )
            plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color = txt_color ) 
            cb.ax.yaxis.set_tick_params(color = txt_color)
            cb.ax.xaxis.set_tick_params(color = txt_color)
            cb.ax.xaxis.set_major_formatter(ScalarFormatter())
            cb.ax.xaxis.set_minor_formatter(ScalarFormatter())
            for ax in [ax1]:
                ax.set_aspect("equal")    
            # set colorbar edgecolor 
            cb.outline.set_edgecolor("white")
        ax_xDist.hist(gaussian_field.flatten(),
                    log     =   True ,
                    color   =   "red", 
                    bins    =   50   , 
                    density =   False )
        
        ax_yDist.loglog(kflat, pkflat, marker = ".",ms=0.1,lw=0., color = "blue", alpha = 0.5)
        ax_yDist.loglog(kmag, pk, lw = 2, color = "black", label=r"Data")
        pkan  = pk_function(kmag, kmag[0])#(kmag/kmag[0])**alpha*pk.max()
        
        ax_yDist.loglog(kmag, pkan, lw = 2, color = "red", label="Analyic")#, label = r"$P(k)\propto%.2f$"%alpha)
    
        ax_yDist.set_xlabel(r"$k$")
        ax_yDist.set_ylabel(r"$P(k)/P_0$")
        ax_xDist.set_xlabel(r"$PDF$")
        ax_xDist.set_ylabel(r"$\delta$")
        fig.tight_layout
        ax_yDist.legend()
        fig.savefig(path+"rand_%dD_%d"%(ndim,self.seed))
        
        # Verify the mean and RMS
        self.logger.info("//////////////////////")
        self.logger.info("----------------------")
        self.logger.info("Generated Field Stats:")
        self.logger.info("----------------------")
        self.logger.info(f"Seed: {self.seed:d}")
        self.logger.info(f"Size: {str(size):s}")
        #self.logger.info(f"alpha: {alpha:f}")
        self.logger.info(f"Mean: {np.mean(field):.3e}")
        self.logger.info(f"RMS:  {np.std (field):.3e}")
        self.logger.info("//////////////////////")
        
        with open(path+"file_header_%dD_%05d.dat"%(ndim, self.seed), mode="w") as f:
            f.write(f"integer: i, ")
            f.write(f"double : d, ")
            f.write(f"single : f ")
            f.write("\n")
            f.write(f"Mean: {np.mean(field):.3e}")
            f.write("\n")
            f.write(f"RMS:  {np.std (field):.3e}")
            #f.write("\n")
            #f.write(f"slope:{alpha:.3e}")
            f.write("\n")
            f.write("ndim, i")
            f.write("\n")
            infon = ""
            for i in ["x","y","z"]: infon += "n%s "%i
            f.write("%s, i"%infon)
            f.write("\n")
            for line in field:	
                f.write("First line shape: %s"%str(line.shape))
                break
        return [fig, fig.axes, fig22, fig22.axes] if ndim==3 else [fig, fig.axes]
    

   
   
   
   
if __name__=="__main__":     
        
    ndim = 3
    L    = 150.0
    N    = [128]*3
    zzz = []
    yyy = [] 
    xxx = []
    for i in range(1,int(2)):
        logger.info("computing %d"%i)
        distr = RandomDistribution(
            Boxlen = L,
            seed   = 124,#np.random.randint(1,1e6),
            ndim   = ndim,
            size   = N ,    #[N]*ndim
            
            #lmin = 30,
            #lmax = L*0.5
            )  
            
        def func(k): return k**(-5./3.)
        
        gaussian_field = distr.distribution(pk_function=func, kcorr  = 1./20, sigmak=0.1)
        distr.plot(gaussian_field, path="randDistr/plots", single_cbar= True)
        distr.save(gaussian_field, path="randDistr")
        zzz += [(gaussian_field[ :       , :       , N[2]//2 ].mean())]
        yyy += [(gaussian_field[ :       , N[2]//2 , :       ].mean())]
        xxx += [(gaussian_field[ N[2]//2 , :       , :       ].mean())]
    plttt=plt.figure(10394)
    axxx=plttt.add_subplot(111)
    axxx.hist(zzz, bins=20, label="z", alpha =0.5)
    axxx.hist(yyy, bins=20, label="y", alpha =0.5)
    axxx.hist(xxx, bins=20, label="x", alpha =0.5)
    axxx.legend()
    plttt.savefig("ffhdf")