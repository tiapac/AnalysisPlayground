from utils.myparser import myparser
from utils.utils import *
import numpy as np
import os
from scipy.io import FortranFile
import time as Time
import constants as cc
from paper_settings import PaperSettings
import warnings
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from utils.geometry import coordinates
from utils.stat import weighted_profile
from FastPlots import dataset
from utils.ErrorHandling import Errors
from utils.logger import setup_logger

warnings.filterwarnings("ignore")
FortranFile.read_vector = FortranFile.read_record
nonolist =  [213, 341, 74]

def compute_2d_ellipse_parameters(xx_flat, yy_flat, m, cond):
    """
    Computes the center of mass, inertia tensor, and equivalent ellipse parameters
    for a 2D mass distribution.
    
    Parameters:
      datas: dict with keys "mass", "xx", "yy". The arrays must all have the same shape.
      cond:  boolean array (of the same flattened length) indicating which points to include.
      
    Returns:
      xm, ym      : Center of mass coordinates.
      a, b        : Semi-axis lengths of the equivalent ellipse. (Note: a <= b,
                    so here a is the semi-minor and b the semi-major axis.)
      theta       : Orientation (in radians) of the major axis measured from the x-axis.
      I           : 2x2 inertia tensor about the center of mass.
      eigenvalues : Eigenvalues of the inertia tensor.
      eigenvectors: Corresponding eigenvectors (columns); the eigenvector of b is used for theta.
    """
    # Flatten the data arrays.

    
    # Apply the condition to select points.
    m_flat = m[cond]
    xf = xx_flat[cond]
    yf = yy_flat[cond]
    
    # Print lengths for debugging.
    print("Number of mass elements:", len(m_flat), "and number of coordinate elements:", len(xf))
    
    # Compute the total mass.
    M = np.sum(m_flat)
    
    # Compute center of mass.
    xm = np.sum(m_flat * xf) / M
    ym = np.sum(m_flat * yf) / M
    
    # Compute the coordinate differences relative to the center of mass.
    dx = xf - xm
    dy = yf - ym
    
    # Compute the second central moments. For a 2D lamina, one common definition of the
    # inertia tensor about the center of mass is:
    #   I_xx = Σ m * (dy²)
    #   I_yy = Σ m * (dx²)
    #   I_xy = -Σ m * (dx*dy)
    I_xx = np.sum(m_flat * dy**2)
    I_yy = np.sum(m_flat * dx**2)
    I_xy = -np.sum(m_flat * dx * dy)
    
    # Form the 2x2 inertia tensor.
    I = np.array([[I_xx, I_xy],
                  [I_xy, I_yy]])
    
    # Diagonalize the inertia tensor to find the principal moments.
    eigenvalues, eigenvectors = np.linalg.eigh(I)
    
    # Sort the eigenvalues and eigenvectors so that eigenvalues[0] is the smaller one.
    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    
    # For a uniform elliptical lamina, the relation between the principal moment I_i and
    # the corresponding semi-axis (say, axis_i) is:
    #       I_i = (M/4) * (axis_i)²
    # Hence, we can recover:
    #       axis_i = sqrt(4 * I_i / M)
    # (Here, by convention we take a <= b so that a corresponds to the smaller eigenvalue.)
    a = np.sqrt(4 * eigenvalues[0] / M)
    b = np.sqrt(4 * eigenvalues[1] / M)
    
    # Determine the orientation angle.
    # We use the eigenvector corresponding to the larger eigenvalue (b) for the major axis.
    # The angle theta (in radians) measured from the x-axis is given by:
    theta = np.arctan2(eigenvectors[1, 1], eigenvectors[0, 1])
    
    return xm, ym, a, b, theta, I, eigenvalues, eigenvectors



def makeplot(self, num):
        if num in nonolist: return
        for ind, tmpkind in enumerate(self.kind):
            
            if self.plot: ax = self.axs[0] if self.args.oneplot else self.axs[ind]
            
            tmpkind = self.kind[ind]
            tmpprojection = self.projection[ind]
            loggx = self.logx[ind]
            loggy = self.logy[ind]
            tmpagainst = self.against[ind]
            tmed = self.args.tquant[0]
            multiindex = 0 if self.multitime is None else self.multitime.index(num)
            current_axis = self.axes[tmpprojection]
            pswitch = tmpkind == "press" and self.args.pswitch

            datas = dataset(num=num, proj=tmpprojection, path=self.path)
            data = datas[tmpkind]
            if self.axrotate[current_axis]:
                data = np.rot90(data)
            time = datas.time * datas.cd.scale_t / cc.kyr

            if tmpkind == "pmag" or pswitch:
                data = data / cc.kB
            if pswitch:
                data2 = datas["pmag"] / cc.kB
            if self.args.normalize is not None:
                data /= abs(data.max()) if self.args.normalize == "max" else float(self.args.normalize)

            data_flat = data.flatten()
            xx_flat = datas["x"].flatten()
            yy_flat = datas["y"].flatten()
            coord = coordinates(xx_flat, yy_flat)
            radius_flat = datas[tmpagainst].flatten()

            if self.args.weight is None:
                mass_flat = np.ones_like(data_flat)
            else:
                mass_flat = datas[self.args.weight].flatten()

            treshold = 2
            vmin = self.clim[ind][0] if self.clim[ind][0] is not None else data.min() * 1.0
            vmax = self.clim[ind][1] if self.clim[ind][1] is not None else data.max() * 0.9
            myaxes = find_axes(self.axes[tmpprojection])

            match tmpprojection:
                case 1:
                    self.logger.info("just a projection. Ignore for now.")
                    continue
                case 2 | 3 | 4:
                    if pswitch:
                        data_flat2P = data2.flatten()

                    if self.plot and self.multitime is None:
                        ax.clear()

                    cond = secure_eval(self.args.cond, {}, {'datas': datas,
                                                            "cc":cc,
                                                            "flat": lambda x: (datas[x] if not self.axrotate[current_axis] else np.rot90(datas[x])).flatten()})
                    # print(len(datas["pmag"].flatten()[datas["pmag"].flatten()>1e3]))
                    # quit()
                    cond1 = coord.filter_coordinates(np.pi / 5.)
                    cond2 = coord.filter_coordinates(np.pi / 5., inverted=True)

                    qmed1, cell_centers = weighted_profile(data_flat[cond1 & cond], radius_flat[cond1 & cond], mass_flat[cond1 & cond], bins=self.args.bins)
                    qmed2, cell_centers = weighted_profile(data_flat[cond2 & cond], radius_flat[cond2 & cond], mass_flat[cond2 & cond], bins=self.args.bins)
                    qmed, cell_centers = weighted_profile(data_flat[cond], radius_flat[cond], mass_flat[cond], bins=self.args.bins)

                    if self.args.normalize == "max":
                        qmed1 /= qmed1.max()
                        qmed2 /= qmed2.max()
                        qmed /= qmed.max()

                    if self.args.normalize is not None:
                        qmed /= qmed.max() if self.args.normalize == "max" else float(self.args.normalize)
                        qmed1 /= qmed1.max() if self.args.normalize == "max" else float(self.args.normalize)
                        qmed2 /= qmed2.max() if self.args.normalize == "max" else float(self.args.normalize)

                    if self.analyse or self.args.vline:
                        if not self.args.inertia:
                            pos_med = find_shock_reverse(qmed, cell_centers, treshold, tmed, invert=self.args.invert, smaller=self.args.smaller)
                            a = find_shock_reverse(qmed1, cell_centers, treshold, tmed, invert=self.args.invert, smaller=self.args.smaller)
                            b = find_shock_reverse(qmed2, cell_centers, treshold, tmed, invert=self.args.invert, smaller=self.args.smaller)
                            xm, ym = 0, 0
                        else:
                            xm, ym, b, a, theta, I, eigenvalues, eigenvectors = compute_2d_ellipse_parameters(xx_flat, yy_flat, mass_flat, cond)
                            pos_med = (a+b)*0.5
                    CommArgs = dict(lw=1)
                    if self.plot and not self.args.onlysph:
                        ycolor0 = "crimson" if self.multitime is None else self.linecmap[multiindex]
                        ycolor1 = "midnightblue" if self.multitime is None else self.linecmap[multiindex]
                        label = myaxes[0] if self.multitime is None else rf"{time:.0f} kyr"
                        ax.plot(cell_centers, qmed1, label=label, color=ycolor0, **CommArgs, linestyle="solid")

                        label = myaxes[1] if self.multitime is None else None
                        ax.plot(cell_centers, qmed2, label=label, color=ycolor1, **CommArgs, linestyle="dashed")

                        if pswitch:
                            qmed1, cell_centers = weighted_profile(data_flat2P[cond1 & cond], radius_flat[cond1], mass_flat[cond1], bins=self.args.bins)
                            qmed2, cell_centers = weighted_profile(data_flat2P[cond2 & cond], radius_flat[cond2], mass_flat[cond2], bins=self.args.bins)
                            ax.plot(cell_centers, qmed1, label=myaxes[0], color="crimson", linestyle="dashed", **CommArgs)
                            ax.plot(cell_centers, qmed2, label=myaxes[1], color="midnightblue", linestyle="dashed", **CommArgs)
                        if self.args.vline:
                            ax.axvline(a, color="black")
                            ax.axvline(b, color="black")

                    if self.plot and not self.args.nosph:
                        ycolor = "black" if self.multitime is None else self.linecmap[multiindex]
                        label = "spherical" if self.multitime is None else rf"{time:.0f} kyr"
                        if self.args.fit:
                            label = ""
                        if not self.args.onlyfit:
                            ax.plot(cell_centers, qmed, label=label, color=ycolor, **CommArgs)
                        coeffs = np.polyfit(cell_centers, qmed, 1)
                        linear_fit = np.poly1d(coeffs)
                        if self.args.fit:
                            ax.plot(cell_centers, linear_fit(cell_centers), label=f"m, q = {coeffs[0]:.2f}, {coeffs[1]:.2f} {time:.0f} kyr",
                                    color=ycolor, **CommArgs)
                        if self.args.vline:
                            ax.axvline(pos_med, color="black")
                        if pswitch:
                            qmedc, cell_centersc = weighted_profile(data_flat2P, radius_flat, mass_flat, bins=self.args.bins)
                            ax.plot(cell_centersc, qmedc, label=label, color=ycolor, linestyle="dashed", **CommArgs)

                    self.logger.info(f"Radius done for {num:05d}, time {time:.3f} on proj {tmpprojection:d} for {tmpkind:s}, axis {self.axes[tmpprojection]:s}-")

                    if self.analyse or self.args.vline:
                        self.logger.info(f"Radii found {pos_med:.4f}, {min(a, b):.4f}, {max(a, b):.4f}")
                        if self.analyse:
                            if self.args.fit:
                                ordered = (num, time, pos_med, min(a, b), max(a, b), coeffs[0], coeffs[1])
                            else:
                                ordered = (num, time, pos_med, min(a, b), max(a, b))
                            with open(self.files[current_axis], mode="a") as myfile:
                                myfile.write(("%d " + "%f " * (len(ordered) - 1) + "\n") % ordered)

                    if self.plot:
                        
                        #plt.setp(ax.get_yticklabels()[0], visible=False)  
                        
                        
                        if not self.args.polar:
                            ax.set_xlim(cell_centers.min(), cell_centers.max())
                        
                  
                        ax.set_xscale("log" if loggx else "linear")
                        ax.set_yscale("log" if loggy else "linear")
                        if self.args.normalize is None:
                            label = fr"{str(self.myunits.symbols[tmpkind])}" + str(self.myunits.units[tmpkind])
                        else:
                            label = fr"{str(self.myunits.symbols[tmpkind])}/{str(self.myunits.symbols[tmpkind])}" + r"$_\mathrm{max}$"
                        ax.set_ylabel(label, **self.paper.figfontstyle)
                        if tmpkind == "plasma_beta" and self.multitime:
                            colors = plt.cm.viridis(np.linspace(0, 1, len(self.multitime)))
                            ycolor = colors[multiindex]

                            if not self.args.nosph:
                                self.axins.plot(cell_centers, qmed, color=ycolor)

                            if not self.args.onlysph:
                                self.axins.plot(cell_centers, qmed1, color=ycolor0, linestyle = "solid")
                                self.axins.plot(cell_centers, qmed2, color=ycolor1, linestyle = "dashed")
                            self.axins.set_yscale("log")
                            ax.indicate_inset_zoom(self.axins, edgecolor="black")

                        if not self.multitime:
                            ax.set_xlabel(fr"{str(self.myunits.symbols[tmpagainst])}" +
                                          str(self.myunits.units[tmpagainst]), **self.paper.figfontstyle)
                            if not self.args.nolabel:
                                ax.legend()
                        else:
                            if not self.args.nolabel:
                                # if ind == 0 and num == self.multitime[-1]:
                                #     ax.legend()
                                if ind == 2:
                                    ax.legend()
                            if ind < self.nplots - 1:
                                ax.set_xticks([])
                            
                            else:
                                ax.set_xlabel(fr"{str(self.myunits.symbols[tmpagainst])}" +
                                              str(self.myunits.units[tmpagainst]), **self.paper.figfontstyle)
                        
                        ax.set_yticks(ax.get_yticks()[:-2])

                        if self.args.xlim and not self.args.polar:
                            ax.set_xlim(self.args.xlim[0], self.args.xlim[1])
                        ax.set_ylim(vmin, vmax)
                        if self.args.sideplot:
                            toplot = data if self.args.splot is None else datas[self.args.splot]
                            splot= tmpkind if self.args.splot is None else self.args.splot
                            myaxs1 = self.axs1[ind][multiindex]
                            myaxs1.clear()
                            myaxs1.imshow(toplot, extent=(xx_flat.min(), xx_flat.max(), yy_flat.min(), yy_flat.max()),
                                          cmap=self.paper.cmaps.colormaps[splot],
                                          norm=LogNorm())
                            if self.args.vline:
                                t = np.linspace(0, 2 * np.pi, 100)
                                myaxs1.plot(xm + a * np.cos(t), ym + b * np.sin(t), color ="black")
                                if not self.args.nosph:  myaxs1.plot(xm + pos_med * np.cos(t), ym + pos_med * np.sin(t), linestyle="dashed")
                            myaxs1.set_xticks([])
                            myaxs1.set_yticks([])
                        

        if self.plot:
            os.makedirs(self.savepath, exist_ok=True)
            name = f"{self.savepath}/profile_{num:05d}"
            self.fig.tight_layout( rect=[0, 0.03, 1, 0.95])
            self.fig.subplots_adjust( hspace = 0.1, wspace = 0.1 )
            if not self.multitime:
                self.fig.savefig(name)
            if self.args.show:
                plt.show()
class FastProfile: 
    def __init__(self, args):
        self.args = args
        self.logger = setup_logger("FastProfile")
        self.paper = PaperSettings()
        self.myunits = self.paper.units
        self.ave = args.noave
        self.path = args.path
        self.plot = args.plot
        self.start = args.n1
        self.end = args.n2
        self.kind = args.quant
        self.projection = args.proj
        self.against = args.against
        self.multitime = self.get_multitime()
        self.nplots = len(self.kind)
        self.savepath = self.get_savepath()
        self.fsavepath = self.get_fsavepath()
        self.axrotate = dict(x=args.rotatex, y=args.rotatey, z=args.rotatez)
        self.cmaps = [self.paper.cmaps.colormaps[k] for k in self.kind]
        self.logx = [not args.nologx] * self.nplots
        self.logy = [not args.nology] * self.nplots
        self.clim = self.get_clim()
        self.initime = Time.time()
        self.labels = [None for lab in range(len(self.kind)) ]
        self.sideplot = args.sideplot
        if self.multitime is not None: self.linecmap = plt.cm.viridis(np.linspace(0,1,self.multilen))

        self.fig, self.axs, self.axins, self.axs1 = self.create_figure() #if args.plot else (None, None, None, None)
        
        self.analyse = args.analyse
        self.namebase = "axes_length_"
        
        self.xyz = "xyz"
        self.axes_string = "xzyx"
        self.axes    = { pj:self.axes_string[pj-1] for pj in self.projection}
        self.files = self.get_files()
    def get_multitime(self):
        if self.args.multitime is not None:
            if self.args.multitime[0] == "range":
                return list(np.arange(int(self.args.multitime[1]),
                                      int(self.args.multitime[2]),
                                      int(self.args.multitime[3])))
            else:
                return sorted([int(n) for n in self.args.multitime])
        return None

    def get_savepath(self):
        basedirname = "profiles/"
        axes_string = "xzyx"
        axes = {pj: axes_string[pj-1] for pj in self.projection}
        self.jointedNames = jointedNames = "".join([f"{k}_" if idx < self.nplots-1 else f"{k}" for idx, k in enumerate(self.kind)])
        self.jointedAgNames = jointedAgNames = "".join([f"{k}_" if idx < self.nplots-1 else f"{k}" for idx, k in enumerate(self.against)])

        if self.multitime is  None:
            savepath = self.args.dir+f"/{basedirname}{jointedAgNames}/{jointedNames}/los_{"".join([str(axes[i]) for i in self.projection])}/"

        elif self.multitime is not None:
            self.multilen = multilen  = len(self.multitime)
            if "range" in self.multitime:
                self.jointedOuts = "".join([f"{k}_" if idx < multilen-1 else f"{k}" for idx,k in enumerate(self.multitime)])
            else:    
                self.jointedOuts = "".join([f"{k:04d}_" if idx < multilen-1 else f"{k}" for idx,k in enumerate(self.multitime)])
            self.jointedLos = "".join([str(axes[i]) for i in self.projection])
            savepath = self.args.dir+f"/{basedirname}Multitime/"
        return savepath


    def get_fsavepath(self):
        basedirname = "profiles/"
        if self.args.fsavepath is None:
            return self.args.dir + f"/dataAxis/"
        else:
            return basedirname + self.args.fsavepath + "/"

    def get_clim(self):
        clim = [[None, None] for i in range(self.nplots)]

        if self.args.clims is not None:
            for i, g in enumerate(self.args.clims):
                for j in [0, 1]: clim[i][j] = None if( g[j] == "auto" or g[j] == "a") else float(g[j])
        return clim

    def checktime(self, msg=""):
        self.logger.info("--- %s %.2e seconds ---" % (msg, (Time.time() - self.initime)))

    def create_figure(self):
        if self.plot:
            import matplotlib.pyplot as plt
            from matplotlib.colors import LogNorm
            from matplotlib import collections  as mc
            plt.rcParams.update({"text.usetex": True, **self.paper.mplfigfontstyle})
            nrows = self.nplots
            if self.args.oneplot: nrows = 1
            nsideplots = 1 if self.multitime is None else self.multilen
            ncols = nsideplots+1 if self.args.sideplot else 1

            fig, _axs = plt.subplots(nrows = nrows, 
            ncols = ncols,
            subplot_kw={'projection': 'polar'} if self.args.polar else None,
            **{"figsize":(self.paper.onecol * ncols, self.paper.onecol * nrows)})
            if ncols ==1 and nrows ==1:
                axs = [_axs]
            elif ncols ==2 and nrows==1:
                axs = [_axs[0]]
                if self.sideplot: axs1 = [[_axs[1]]]
                
                
            else: 

                if nrows ==1: 
                    axs = [_axs[0]] if self.sideplot else _axs.flatten() 
                    if self.sideplot: 
                        if self.multitime is None: 
                            axs1 = [_axs[1]]
                        else:
                            
                            axs1 = [_axs[1:self.multilen+1 ]]
                else: 
                    axs = _axs[:,0] if self.sideplot else _axs.flatten()
                    if self.sideplot: 
                        if self.multitime is None: 
                            axs1 = [_axs[:,1]]
                        else:
                            
                            axs1 = [_axs[:,1:self.multilen+1 ]]
        else: 
            axs = range(self.nplots)

        if "plasma_beta" in self.kind:
            x1, x2, y1, y2 = 30, 50, 5e-2, 5  # subregion of the original image
            axins = axs[self.nplots-1].inset_axes([0.5, 0.5, 0.47, 0.47],xlim=(x1, x2), ylim=(y1, y2))#, xticklabels=[], yticklabels=[])
            first = False
        else: 
            axins = None
        if not self.sideplot: axs1 = None
        return fig, axs, axins, axs1

    def get_files(self):
        if self.analyse:
            filey = self.fsavepath + self.namebase + "y.dat"
            filex = self.fsavepath + self.namebase + "x.dat"
            filez = self.fsavepath + self.namebase + "z.dat"
            tfiles = dict(x=filex, y=filey, z=filez)
            files = dict()
            for pj in self.projection:
                current_axis = self.axes[pj]
                files.update({current_axis: tfiles[current_axis]})
            os.makedirs(self.fsavepath, exist_ok=True)
            if self.args.clean:
                for file in files.values():
                    if os.path.exists(file):
                        os.remove(file)
            return files
        return None

    

    def run(self):
        if not self.multitime:
            snaplist = range(self.start, self.end, self.args.step)
            if self.args.cpus is not None:
                from itertools import repeat
                from multiprocessing import Pool
                with Pool(processes=self.args.cpus) as p:
                    results = p.starmap_async(makeplot, zip(repeat(self),snaplist))
                    results.get()
                    p.close()
            else:
                for num in range(self.start, self.end, self.args.step):
                    makeplot(self,num)
        else:
            if self.args.cpus is not None:
                raise Errors.IllegalArgumentError("Cant use multiprocessing with -multitime.")
            name = self.savepath + f"{self.jointedNames}_times_{self.jointedOuts}_los_{self.jointedLos}"
            for num in list(self.multitime):
                makeplot(self,num)

            if self.plot:
                psize = self.paper.psize
                self.fig.set_size_inches(*psize[self.args.impage])
                self.fig.tight_layout( rect=[0, 0.03, 1, 0.95])
                #self.fig.subplots_adjust( hspace = 0.0, wspace = 0.0 )
                self.fig.savefig(name, dpi = self.paper.figure.dpi)


def MakeFastProfile(parser = myparser(prog = "Fastplosts"), subc_args = None):
    import numpy as np
    import os
    #from cython_fortran_file import FortranFile
    from scipy.io import FortranFile
    FortranFile.read_vector = FortranFile.read_record
    import time as Time
    import constants as cc

    from paper_settings import PaperSettings
    import warnings
    warnings.filterwarnings("ignore")
    import matplotlib.pyplot as plt
    # from matplotlib.ticker import ScalarFormatter, LogFormatter
    from matplotlib.colors import LogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # from mpl_toolkits.axes_grid1 import make_axes_locatable
    # from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

    from utils.geometry import coordinates
    from utils.stat import weighted_profile
    from FastPlots import dataset
    from utils.ErrorHandling import Errors
    from utils.logger import setup_logger


    logger = setup_logger("FastProfile")
    #sfilter_coordinates = coordinates.filter_coordinates

    paper = PaperSettings()
    myunits = paper.units


    plt.rcParams.update({
                    "text.usetex": True, **paper.mplfigfontstyle
                    #"font.family": "Helvetica"
                    })
    # parser = myparser(prog = "FastProfiles", overwrite=True)
    
    parser.add_argument("-n1", type=int, default=None, help="Enter at least one output number")

    parser.add_argument("-q", "--quant", required=True, action="append", type=str, help="Specify the quantity to be plotted")
    parser.add_argument("-ag", "--against", nargs="+", default=["radius"], type=str, help="Specify the variable to plot against")
    parser.add_argument("-fs", "--fsavepath", default=None, type=str, help="Specify the save path for the output files")
    parser.add_argument("-s", "--step", type=int, default=1, help="Specify the step size for the output numbers")

    parser.add_argument("-p", "--proj", type=int, required=True, action="append", help="Specify the projection to be used")
    parser.add_argument("-b", "--bins", type=int, required=False, default=50, help="Specify the number of bins for the profile")
    parser.add_argument("-sideplot", action="store_true", default=False, help="Enable side plot")
    parser.add_argument("-inertia", action="store_true", default=False, help="Enable side plot")
    parser.add_argument("-noave", action="store_false", default=True, help="Disable averaging")
    parser.add_argument("-analyse", action="store_true", default=False, help="Enable analysis mode")
    parser.add_argument("-plot", action="store_true", default=False, help="Enable plotting")
    parser.add_argument("-onlysph", action="store_true", default=False, help="Plot only spherical data")
    parser.add_argument("-nosph", action="store_true", default=False, help="Plot only spherical data")
    parser.add_argument("-threeway", action="store_true", default=False, help="Enable three-way plotting")
    parser.add_argument("-path", default="./", help="Specify the path to the data files")
    parser.add_argument("-multitime", default=None, help="List of snapshots", nargs="+")
    parser.add_argument("-impage", default=1, type=int, help="Specify the image page size")
    parser.add_argument("-name", default=None, help="Specify the name for the output files")
    parser.add_argument("-xlim", default=None, nargs=2, type=float, help="Specify the x-axis limits")
    parser.add_argument("-normalize", default=None, help="Specify normalization method")
    parser.add_argument("-pswitch", action="store_true", default=False, help="Enable pressure switch")
    parser.add_argument("-vline", action="store_true", default=False, help="Enable vertical line in plots")
    parser.add_argument("-clean", action="store_true", default=False, help="Clean existing files")
    parser.add_argument("-invert", action="store_true", default=False, help="Invert the profile")
    parser.add_argument("-tquant", nargs=1, default=[2.0e4], type=float, help="Specify the temperature quantile")
    parser.add_argument("-c", "--clims", action="append", nargs=2, default=None, help="Specify color limits")
    parser.add_argument("-nlx", "--nologx", action="store_true", default=False, help="Disable log scale for x-axis")
    parser.add_argument("-nly", "--nology", action="store_true", default=False, help="Disable log scale for y-axis")
    parser.add_argument("-nl", "--nolog", action="store_true", default=False, help="Disable log scale for both axes")
    parser.add_argument("-rx", "--rotatex", action="store_true", default=False, help="Disable log scale for both axes")
    parser.add_argument("-ry", "--rotatey", action="store_true", default=False, help="Disable log scale for both axes")
    parser.add_argument("-rz", "--rotatez", action="store_true", default=False, help="Disable log scale for both axes")
    parser.add_argument("-w", "--weight", default=None, help="Disable weighting")
    parser.add_argument("-splot", default=None, help="Disable weighting")
    parser.add_argument("-nolabel", action="store_true", default=False, help="Disable labels")
    parser.add_argument("-pl", "--polar", action="store_true", default=False, help="Enable polar projection")
    parser.add_argument("-fit", action="store_true", default=False, help="Enable linear fit")
    parser.add_argument("-onlyfit", action="store_true", default=False, help="Plot only the linear fit")
    parser.add_argument("-show", action="store_true", default=False, help="Plot only the linear fit")
    parser.add_argument("-oneplot", action="store_true", default=False, help="Plot only the linear fit")
    parser.add_argument("-smaller", action="store_true", default=False, help="Plot only the linear fit")
    parser.add_argument("-cond", type=str, default="True", help="Specify the condition for filtering data")
    parser.add_default()
    args = parser.parse_args(subc_args,fix_n1n2 = False, fixdir=True)
    if args.n1 is not None: 
        args.n1,args.n2 = parser._fix_n1n2(args)

    if args.nolog: 
        args.nologx = True
        args.nology = True
    if args.n1 is None and args.multitime is None: 
        raise Errors.IllegalArgumentError("-n1 must be specifield if -multitime is not.") 


    kind       = args.quant          #["press","press","press"]
    projection = args.proj      #[4, 4, 2, 2]
    against    = args.against 
    if len(projection)<len(kind):
        for i in range(len(projection)-1,len(kind)-len(projection)):
            projection.append(projection[i-1])
            
    if len(kind)<len(projection):
        for i in range(len(kind)-1,len(projection)-len(kind)):
            kind.append(kind[i-1])

    if len(against)<len(projection):
        for i in range(len(against)-1,len(projection)-len(against)):
            against.append(against[i-1])
    args.quant   = kind                  #["press","press","press"]
    args.proj    = projection        #[4, 4, 2, 2]
    args.against = against      


    if args.threeway:
        assert len(kind)==1, "threeway only allowed with one kind"
        args.projection = [2,3,4] 
        
    dd = FastProfile(args)
    # dd.clim   = [
    #             [2e-3 , 2e1 ], [2e-3, 2e1  ],
    #             [2e2  , 1e8 ], [2e2 , 1e8  ],
    #             [1e-9 , 1e-4], [1e-9, 1e-4 ],
    #             [-2   , 5   ], [-2  , 5    ]
    #             ]
    dd.run()