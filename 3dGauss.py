import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import paper_settings as pg 
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter, NullFormatter
import os
#from pyevtk.hl import gridToVTK

#from scipy.io import FortranFile
from cython_fortran_file import FortranFile
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
paper = pg.PaperSettings()

plt.rcParams.update({
    "text.usetex": True, **paper.mplfigfontstyle
    #"font.family": "Helvetica"
    })
themes = {"div" : dict(cmap = "RdBu", color="black"),
          "dark": dict(cmap = "cividis", color="white")}
theme    = themes["div"]
txt_color = theme["color"]
cmap      = theme["cmap"]





def PRN(size: tuple|list,
        mean : float|int,
        rms: float|int  ,
        correlation_length: float|int,
        dx : float|int = 1 ,
        alpha : float|int = 5./3.,
        seed  : float|int = 1,
        kmin:float=None,
        kmax:float=None,
        bins = 100 
        ):
    
    np.random.seed(seed)

    # computie dimensionality 
    ndim   = len(size)    
    
    # create the frequencies from the signal's number with np.fttfreq 

    # this is equivalent to     ki = [kx = np.fft.fftfreq(nx),...,   kn = np.fft.fftfreq(nn)]
    ki = [  np.fft.fftfreq(n, d = dx).reshape([-1 if i == j else 1 for j in range(ndim)]) for i, n in enumerate(size) ]
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
    if kmin is None and kmax is None:
        pk = k**alpha
        
        # now we normalize the specturm base on the grid volume     
        pk /= np.sum(pk)     # over the sum 
        pk *= np.prod(size)  # and scale by total number of grid points

        # Now that we have the specturm, we now want to randomize the complex amplutudes of the signal. 
        # therefore, we generate a gaussian noise
            
        noise = np.random.normal(0, 1, size) + 1j * np.random.normal(0, 1, size)
         # we now randomize our nodes and thake the amplitudes, nowing that
        # |A(k)|**2 = P(k) -> |A(k)| = (P(k))**1/2 
        fourier_field = noise * np.sqrt(pk)
    else:
    ### do computations on a specific log-spaced range of k
    # mmmh....
    
    
        kbins = np.logspace(np.log10(kmin), np.log10(kmax), bins)
        pkbins = kbins**alpha

        # Assign bin indices to k values
        k_indices = np.digitize(k.flat, kbins, right=True)

        # Initialize power spectrum
        pk = np.zeros_like(k)

        # Assign power spectrum to bins
        for i, pkbin in enumerate(pkbins):
            pk.flat[k_indices == i] = pkbin

        # Normalize the power spectrum
        pk /= np.sum(pk)
        pk *= np.prod(size)

        # Generate complex Gaussian noise
        noise = np.random.normal(0, 1, size) + 1j * np.random.normal(0, 1, size)
        #noise[k>kmax]=0.0
        # Apply the power spectrum
        fourier_field = noise * np.sqrt(pk)
        
    
    
        
   
    
    
    
    
    
    
    # back to real space with the inverse of the fourier distribution
    # but we only keep the real component 
    field = np.fft.ifftn(fourier_field).real  
    
    # normalize again the field 
    field = field / np.std(field)
    
    # rescale the field 
    field *= rms
    
    # shift to the desired mean 
    field += mean
    return field


def compute_pk(field, kmax = None, kmin = None ):
    
    size = field.shape
    gridvolume = np.prod(size)
   
    # compute the fourier transform
    fourier_field = np.abs(np.fft.fftn(field))**2  # Take the squared magnitude

    #print(fourier_field.shape)
    ki = [  np.fft.fftfreq(n).reshape([-1 if i == j else 1 for j in range(ndim)]) for i, n in enumerate(size) ]
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
    

# Parameters

save_plot_path = "./plottest/" 
save_data_path = "./data/"
plot = True
L 	  = 150.0
alpha = -5./3.#-5./3.

seed   = 7
level = 8
# dimensionality of the system
ndim   = 3
# resolution 
N      = int(2**level)


lmin = 20#2*np.pi
lmax = L/2
kmin = 2*np.pi/lmax#1./10#/size.min()
kmax = 2*np.pi/lmin

# Size of the field

size  = [N] * ndim  

size = np.array(size, dtype="i")

dx 	  = L / N 
coord = [np.arange(N)*dx - 0.5 * dx]*ndim

mean  = 0.       # Desired mean of the field
rms   = 0.1        # Desired RMS of the field
lcorr = 2
correlation_length = lcorr * N / L

# Generate the field
fignum = 1
import matplotlib.gridspec as gridspec
if not os.path.exists(save_plot_path): os.makedirs(save_plot_path)
if not os.path.exists(save_data_path): os.makedirs(save_data_path)
        
first = True

for seed in range(10,11): 
    gaussian_field = PRN(size      ,
                        mean      , 
                        rms       ,
                        lcorr     ,
                        alpha  = alpha,
                        dx 	= dx,
                        seed   = seed,
                        kmin = kmin,#/size.min()
                        kmax = kmax #size.max()*20
    
                        )
    kmag, pk, kflat, pkflat = compute_pk(gaussian_field, kmin = kmin,kmax = kmax )

    if plot:
        if first:
            fig      = plt.figure(fignum)
            gs       = gridspec.GridSpec(3, 3)
            ax1      = fig.add_subplot(gs[1:3, :2])
            ax_xDist = fig.add_subplot(gs[0,   :2] )
            cbaxes = inset_axes( ax1, width = "70%" , height      = "5%", loc = "upper center")    
            
            fig.set_size_inches(2*paper.onecol,2*paper.onecol)
            
            ax_yDist = fig.add_subplot(gs[1:3, 2])
            first = False
        else: 
            for ax in fig.axes: ax.cla()
        #fig, [ax2,ax3] = plt.subplots(2,1, figsize = (4,24))
        if ndim > 1: 
            #ext = [gaussian_field.min(),gaussian_field.max()] #max(abs(gaussian_field.min()), abs(gaussian_field.max()))
            
            if ndim == 2:
                h = ax1.pcolormesh(coord[0], coord[1], gaussian_field, cmap = cmap)#, vmin=vmin, vmax=vmax)

                #h = ax1.pcolormesh(coord[0], coord[1], gaussian_field, cmap = cmap, vmin=-ext, vmax=ext)

            else:
                
                cut = size[0]//2                
                triplot=True
                toplot = gaussian_field[:, :, cut]
                ext = [toplot.min(),toplot.max()]
                vmin,vmax = ext

                h = ax1.pcolormesh(coord[0], coord[1],toplot , cmap = cmap, vmin=vmin, vmax=vmax)

#                h = ax1.pcolormesh(coord[0], coord[1], gaussian_field[:, :, N//2], cmap = cmap, vmin=-ext, vmax=ext )
                if triplot:
                    fig22,axx = plt.subplots(1,3, figsize=(12, 4))
                    #[ [ax11, ax12], [ax21 , ax22]]  = axx
                    [ ax11, ax12,ax21 ]  = axx
                    #toplot = np.sum(gaussian_field, axis=0)
                    toplot = gaussian_field[:, :, cut]
                    ext = [toplot.min(),toplot.max()]
                    ax11.pcolormesh(coord[0], coord[1],toplot , cmap = cmap, vmin=vmin, vmax=vmax)
                    ax11.set_xlabel("x")
                    ax11.set_ylabel("y")
                    ax11.set_aspect("equal")  
                    
                    #toplot = np.sum(gaussian_field, axis=1)#[:, :, cut]
                    toplot = gaussian_field[:, cut, :]
                    #toplot = gaussian_field[:, :, cut+1]
                    ext = [toplot.min(),toplot.max()]
                    vmin,vmax = ext

                    ax12.pcolormesh(coord[0], coord[2],toplot , cmap = cmap, vmin=vmin, vmax=vmax)
                    ax12.set_xlabel("x")
                    ax12.set_ylabel("z")
                    ax12.set_aspect("equal")  
                    #toplot = np.sum(gaussian_field, axis=2)#[:, :, cut]
                    toplot = gaussian_field[cut, :, :]
                    ext = [toplot.min(),toplot.max()]
                    vmin,vmax = ext
                    #toplot = gaussian_field[:, :, cut+2]
                    
                    
                    ax21.pcolormesh(coord[1], coord[2],toplot , cmap = cmap, vmin=vmin, vmax=vmax)
                    ax21.set_xlabel("y")
                    ax21.set_ylabel("z")
                    ax21.set_aspect("equal")  
                    
                    fig22.tight_layout()
                    fig22.savefig(save_plot_path+"triplot.png")
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
        pkan  = (kmag/kmag[0])**alpha*pk.max()
        
        ax_yDist.loglog(kmag, pkan, lw = 2, color = "red", label = r"$P(k)\propto%.2f$"%alpha)
       
        ax_yDist.set_xlabel(r"$k$")
        ax_yDist.set_ylabel(r"$P(k)/P_0$")
        ax_xDist.set_xlabel(r"$PDF$")
        ax_xDist.set_ylabel(r"$\delta$")
        fig.tight_layout
        ax_yDist.legend()
        fig.savefig(save_plot_path+"rand_%dD_%d"%(ndim,seed))
        # Verify the mean and RMS
        print("//////////////////////")
        print("----------------------")
        print("Generated Field Stats:")
        print("----------------------")
        print(f"Seed: {seed:d}")
        print(f"Size: {str(size):s}")
        print(f"alpha: {alpha:f}")
        print(f"Mean: {np.mean(gaussian_field):.3e}")
        print(f"RMS:  {np.std (gaussian_field):.3e}")
        print("//////////////////////")
        
        with open(save_data_path+"file_header_%dD_%05d.dat"%(ndim, seed), mode="w") as f:
            f.write(f"integer: i, ")
            f.write(f"double : d, ")
            f.write(f"single : f ")
            f.write("\n")
            f.write(f"Mean: {np.mean(gaussian_field):.3e}")
            f.write("\n")
            f.write(f"RMS:  {np.std (gaussian_field):.3e}")
            f.write("\n")
            f.write(f"slope:{alpha:.3e}")
            f.write("\n")
            f.write("ndim, i")
            f.write("\n")
            infon = ""
            for i in ["x","y","z"]: infon += "n%s "%i
            f.write("%s, i"%infon)
            f.write("\n")
            for line in gaussian_field:	
                f.write("First line shape: %s"%str(line.shape))
                break
        #file = FortranFile(save_data_path+"RandomField_%dD_%05d.dat"%(ndim, seed), mode="w")
        file = FortranFile(save_data_path + "RandomField.dat", mode="w")
        
        
        file.write_vector((np.array([ndim]).flatten()))
        file.write_vector(size)

        #for line in np.flip(gaussian_field.T,axis=1):
        #gridToVTK(
        #        "./gaussian_field",
        #        coord[0], coord[1], coord[2],
        #        pointData={"field": gaussian_field}
        #    )
        if ndim==2:
            for line in np.rot90(gaussian_field.T):
               # print(line)
                file.write_vector(line)
        elif ndim==3:
            #for line in gaussian_field.T:#np.rot90(gaussian_field.T):
                #print(line)
                #line = np.rot90(line)
            #    file.write_vector(line)
            mygauss = gaussian_field.copy(order="F")# np.rot90(gaussian_field, axes=(0,2))
            file.write_vector(mygauss.T)
            #for iplane in range(size[2]):
            #    plane = mygauss[:,:, iplane] #np.rot90((mygauss[iplane,: , :]))# "good", without roation
            #    file.write_vector(plane)
#plt.show()