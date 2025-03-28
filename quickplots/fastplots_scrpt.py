from utils.myparser import myparser
import string
alphabet = list(string.ascii_lowercase)
import numpy as np
import matplotlib.pyplot as plt
import os 
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import constants as cc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from paper_settings import PaperSettings
import matplotlib
#from matplotlib.ticker import ScalarFormatter, LogFormatter
from matplotlib.colors import LogNorm
from os.path import dirname
from FastPlots import dataset
# from utils.utils import stopwatch
from utils.logger import setup_logger
from numpy import ma
#from matplotlib import cbook
from matplotlib.colors import Normalize
nonolist =  [213, 341]
def makeplot(self, num, argaxs):
        for ind, ax in enumerate(argaxs):
            ds = dataset(num, proj=self.projection[ind], path=self.path)
            ds.set_loglevel(ds.logger.loglevels.CRITICAL)
            cd = ds.cd
            time = ds.time
            info = ds.info
            data = ds[self.kind[ind]]
            frame_nx, frame_ny = ds.frame_nx, ds.frame_ny
            os.makedirs(dirname(self.savepaths[ind]), exist_ok=True)
            boxlen = info["boxlen"]
            old_info = dict(info)
            update = False
            mask = [old_info[k] != info[k] for k in info.keys()]
            for x, y in zip(mask, info.keys()):
                if x and (y != "nstep_coarse" and y != "time"):
                    update = True
            if self.update_always:
                update = True
            if self.kind[ind] == "pmag" or self.kind[ind] == "press":
                data = data / cc.kB
            if self.kind[ind] == "plasma_beta":
                data = np.log10(data)
                self.myunits.symbols[self.kind[ind]] = r"$\mathrm{Log_{10}}\beta$"
                self.log[ind] = False
            current_axis = self.axes[self.projection[ind]]
            if self.axrotate[current_axis]:
                data = np.rot90(data)
            dataflatten = data.flatten()
            mean = dataflatten.mean()
            self.logger.info("Average data for %s = %.3e, snap %d" % (self.kind[ind], mean, num))
            vmin = self.clim[ind][0] if self.clim[ind][0] is not None else data.min()
            vmax = self.clim[ind][1] if self.clim[ind][1] is not None else data.max()
            self.logger.info("Vmin, vmax %.3e ,%.3e, %.3e ,%.3e, %d" % (vmin, vmax, data.min(), data.max(), num))
            if self.first[ind] or self.args.multitime:
                if not self.args.multitime:
                    ax.clear()
                ax.set_xticks([])
                ax.set_yticks([])
                if self.kind[ind] == "plasma_beta":
                    self.plts[ind] = ax.imshow(data, cmap=self.cmaps[ind], norm=MidPointNorm(0.0))
                else:
                    self.plts[ind] = ax.imshow(data, cmap=self.cmaps[ind], norm=LogNorm(clip=False) if self.log[ind] else None)
                cbargs = dict(width="80%", height="6%", loc="upper center", borderpad=0.5)
                self.plts[ind].set_clim(vmin, vmax)
                cbar2args = dict(orientation="horizontal", extend='both')
                dobar = ind % 2 == 0
                if self.args.multitime is not None:
                    dobar = num == self.args.multitime[0]
                    self.textcolor[self.kind[ind]] = "black"
                    cbargs = dict(width="80%", height="8%", loc="upper center", borderpad=-1.0)
                    cbar2args = dict(location="top", orientation="horizontal", extend='both')
                if dobar:
                    cax = inset_axes(ax, **cbargs)
                    cb = self.fig.colorbar(self.plts[ind], cax=cax, **cbar2args)
                    cb.set_label(fr"{self.myunits.symbols[self.kind[ind]]}" + str(self.myunits.units[self.kind[ind]]),
                                 color=self.textcolor[self.kind[ind]], **self.paper.figfontstyle)
                    cb.ax.yaxis.set_tick_params(color=self.textcolor[self.kind[ind]])
                    cb.ax.xaxis.set_tick_params(colors=self.textcolor[self.kind[ind]], which="both")
                    cb.ax.yaxis.set_tick_params(colors=self.textcolor[self.kind[ind]], which="both")
                    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=self.textcolor[self.kind[ind]])
                    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color=self.textcolor[self.kind[ind]])
                    cb.outline.set_edgecolor(self.textcolor[self.kind[ind]])
                ax.set_xticks([])
                ax.set_yticks([])
                scalebar = AnchoredSizeBar(ax.transData,
                                           self.Lenght_to_Pixels(boxlen / self.paper.figure.barlength[0], boxlen, frame_nx),
                                           r'%.f %s' % (boxlen / self.paper.figure.barlength[0], "pc"),
                                           loc='lower left', pad=0.1, color='white', frameon=False, size_vertical=1)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                if self.args.multitime is not None:
                    if ind == 0:
                        ax.add_artist(scalebar)
                else:
                    if ind == 0:
                        ax.add_artist(scalebar)
                if self.hist[ind]:
                    miniax = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0.5, .05, 0.4, .2),
                                        bbox_transform=ax.transAxes)
                    self.miniplts[ind] = miniax
                    miniax.clear()
                    if self.log[ind]:
                        self.miniplts[ind].set_xscale("log")
                    bins = np.logspace(np.log10(vmin), np.log10(vmax), num=self.nbins[ind])
                    self.miniplts[ind].hist(dataflatten, bins=bins, log=True, color="red")
                    self.miniplts[ind].set_yticks([])
                    self.miniplts[ind].xaxis.set_ticks_position('top')
                    self.miniplts[ind].xaxis.set_label_position('top')
                    miniax.tick_params(axis='both', colors='w', labelcolor=self.textcolor[self.kind[ind]])
                    miniax.spines['bottom'].set_color(self.textcolor[self.kind[ind]])
                    miniax.spines['top'].set_color(self.textcolor[self.kind[ind]])
                    miniax.spines['left'].set_color(self.textcolor[self.kind[ind]])
                    miniax.spines['right'].set_color(self.textcolor[self.kind[ind]])
                    miniax.set_yticks([])
                    miniax.minorticks_off()
                self.first[ind] = False
                # if (len(self.kind) == 1 and self.timestamp) or (self.args.multitime is not None):
                #     if ind == 0:
                #         self.set_timestamp(ax, ind, time, cd, color=self.textcolor[self.kind[ind]])
                # else:
                if self.timestamp and ind == 0:
                    self.set_timestamp(ax, ind, time, cd, color=self.textcolor[self.kind[ind]])
                self.set_panel(ax, ind, color=self.textcolor[self.kind[ind]])
            else:
                for txt in ax.texts:
                    txt.remove()
                if self.timestamp and ind == 0 :
                    self.set_timestamp(ax, ind, time, cd)
                self.set_panel(ax, ind, color=self.textcolor[self.kind[ind]])
                self.plts[ind].set_data(data)
                self.plts[ind].set_clim(vmin, vmax)
                if self.hist[ind]:
                    self.miniplts[ind].clear()
                    self.miniplts[ind].set_xticks([])
                    bins = np.logspace(np.log10(vmin), np.log10(vmax), num=self.nbins[ind])
                    self.miniplts[ind].hist(dataflatten, bins=bins, log=True, color="red")
                    self.miniplts[ind].set_yticks([])
                    self.miniplts[ind].xaxis.set_ticks_position('top')
                    self.miniplts[ind].xaxis.set_label_position('top')
                    if self.log[ind]:
                        self.miniplts[ind].set_xscale("log")
                    self.miniplts[ind].minorticks_off()
        if not self.args.multitime:
            self.logger.info("saving to  " + self.savepaths[ind] + f"/{num:05d}")
            self.fig.tight_layout( rect=[0, 0.03, 1, 0.95])
            self.fig.subplots_adjust( hspace = 0.1, wspace = 0.1 )
            if num not in nonolist: self.fig.savefig(self.savepaths[ind] + f"/{num:05d}", dpi=self.paper.figure.dpi, bbox_inches='tight')
        return argaxs

class MidPointNorm(Normalize):    
        def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
            Normalize.__init__(self,vmin, vmax, clip)
            self.midpoint = midpoint

        def __call__(self, value, clip=None):
            if clip is None:
                clip = self.clip

            result, is_scalar = self.process_value(value)

            self.autoscale_None(result)
            vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

            if not (vmin < midpoint < vmax):
                raise ValueError("midpoint must be between maxvalue and minvalue.")       
            elif vmin == vmax:
                result.fill(0) # Or should it be all masked? Or 0.5?
            elif vmin > vmax:
                raise ValueError("maxvalue must be bigger than minvalue")
            else:
                vmin = float(vmin)
                vmax = float(vmax)
                if clip:
                    mask = ma.getmask(result)
                    result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                    mask=mask)

                # ma division is very slow; we can take a shortcut
                resdat = result.data

                #First scale to -1 to 1 range, than to from 0 to 1.
                resdat -= midpoint            
                resdat[resdat>0] /= abs(vmax - midpoint)            
                resdat[resdat<0] /= abs(vmin - midpoint)

                resdat /= 2.
                resdat += 0.5
                result = ma.array(resdat, mask=result.mask, copy=False)                

            if is_scalar:
                result = result[0]            
            return result

        def inverse(self, value):
            if not self.scaled():
                raise ValueError("Not invertible until scaled")
            vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

            if np.iterable(value):
                val = ma.asarray(value)
                val = 2 * (val-0.5)  
                val[val>0]  *= abs(vmax - midpoint)
                val[val<0] *= abs(vmin - midpoint)
                val += midpoint
                return val
            else:
                val = 2 * (value - 0.5)
                if val < 0: 
                    return  val*abs(vmin-midpoint) + midpoint
                else:
                    return  val*abs(vmax-midpoint) + midpoint
        
class FastPlotter:
    def __init__(self, args):
        print(args)
        self.args = args
        self.logger = setup_logger("FastPlots_script")
        self.logger.setLevel("INFO")
        self.paper = PaperSettings()
        self.myunits = self.paper.units
        self.ave = args.average
        self.fullpage = True
        self.path = args.path
        self.timestamp = args.timestamp
        self.start = args.n1
        self.step = args.step
        self.end = args.n2
        self.kind = args.quant
        self.projection = args.proj
        self.axrotate = dict(x=args.rotatex, y=args.rotatey, z=args.rotatez)
        self.midpoint = False if args.midpoint is None else True
        self.midpoint_v = None if args.midpoint is None else args.midpoint
        self.cmaps = [self.paper.cmaps.colormaps[k] for k in self.kind]
        self.hist = [args.hist] * len(self.kind)
        self.log = [not args.nolog] * len(self.kind)
        self.clim = [[None, None]] * len(self.kind)
        self.axes_string = "xzyx"
        self.axes = {pj: self.axes_string[pj-1] for pj in self.projection}
        self.nbins = [50] * len(self.kind)
        self.textcolor = self.paper.cmaps.textcolor
        self.update_always = True
        self.first = [True] * len(self.kind)
        self.plts = [None] * len(self.kind)
        self.miniplts = [None] * len(self.kind)
        self.savepaths = self.get_savepaths()
        self.fig, self.axs = self.create_figure()
        plt.rcParams.update({"text.usetex": True, **self.paper.mplfigfontstyle}) #"font.family": "Helvetica"})


    def get_savepaths(self):
        nplots = len(self.kind)
        if nplots < 2:
            return [self.args.dir + f"/mymovie_{pj:d}" + f"/{kin}/" for kin, pj in zip(self.kind, self.projection[0:nplots])]
        elif self.args.multitime is not None:
            return [self.args.dir + f"/mymovie_multi_{pj:d}" + f"/{kin}/" for kin, pj in zip(self.kind, self.projection[0:nplots])]
        else:
            jointed = "".join([f"{k}_" if idx < nplots-1 else f"{k}" for idx, k in enumerate(self.kind)])
            return [self.args.dir + f"/mymovie_{"".join([str(i) for i in self.projection])}/{jointed}/"] * nplots

    def create_figure(self):
        nplots = len(self.kind)
        if nplots == 1:
            nrows, ncols = 1, 1
        else:
            nrows = nplots // 2
            ncols = 2 if nplots > 1 else 1
            if nrows + ncols - 1 < nplots:
                ncols += 1
        if self.args.multitime is not None:
            ncols = len(self.args.quant)
            nrows = len(self.args.multitime)
        if self.args.nrows is not None:
            nrows = self.args.nrows
        if self.args.ncols is not None:
            ncols = self.args.ncols
        assert ncols * nrows >= nplots, f"ncols*nrows ({ncols*nrows}) smaller than nplots {nplots}."
        if not self.fullpage:
            fig, iniaxes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(self.paper.onecol * ncols, self.paper.onecol * nrows))
        else:
            fig, iniaxes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(self.paper.twocol * 2, self.paper.twocol * 2))
        if self.args.multitime is not None:
            f = 0.8
            fig, iniaxes = plt.subplots(
                nrows, ncols,
                gridspec_kw=dict(wspace=0.0, hspace=0.0, top=1. - 0.5 / (nrows + 1), bottom=0.5 / (nrows + 1),
                                 left=0.5 / (ncols + 1), right=1 - 0.5 / (ncols + 1)),
                figsize=(self.paper.pagewidth * f, self.paper.pagelength * (f - 0.01)),
                sharey='row', sharex='col'
            )
        axs = [iniaxes] if nplots == 1 else list(iniaxes.flatten())
        fig.set_size_inches(*self.paper.psize[self.args.impage])
        return fig, axs

    def set_timestamp(self, ax, ind, time, cd, color="white"):
        return ax.text(0.01, 0.15, "{:10.0f} kyr".format(time * cd.scale_t / cc.yr / 1e3),
                       fontsize=self.paper.font.body.fontsize, color=color, transform=ax.transAxes)

    def set_panel(self, ax, ind, color="white"):
        return ax.text(0.9, 0.1, f"{alphabet[ind]})", fontsize=self.paper.font.body.fontsize, color=color, transform=ax.transAxes)


    def Lenght_to_Pixels(self, l, boxlen, pixlen):
        return pixlen * l / boxlen

    def run(self):
        snaplist = range(self.start, self.end, self.step)
        if not self.args.multitime:
            if self.args.cpus is not None:
                from itertools import repeat
                from multiprocessing import Pool
                with Pool(processes=self.args.cpus) as p:
                    results = p.starmap_async(makeplot, zip(repeat(self), snaplist, repeat(self.axs)))
                    results.get()
                    p.close()
            else:
                for i in snaplist:
                    makeplot(self, i, self.axs)
        else:
            qnum = len(self.args.quant)
            for j, i in enumerate(self.args.multitime):
                self.axs[j * qnum:(j + 1) * qnum] = makeplot(self, i, self.axs[j * qnum:(j + 1) * qnum])
                self.logger.info("saving to  " + self.savepaths[0] + f"/{i:05d}")
            self.logger.info("saving to  " + self.savepaths[0] + f"/{i:05d}")
            self.fig.tight_layout( rect=[0, 0.03, 1, 0.95])
            self.fig.subplots_adjust( hspace = 0.1, wspace = 0.1 )
            self.fig.savefig(self.savepaths[0] + f"/{i:05d}", self.paper.figure.dpi, bbox_inches='tight')
        quit()




def MakeFastPlot(parser = myparser(prog = "Fastplosts"), subc_args = None):

    parser.add_default()
    parser.add_argument("-q","--quant", type = str, action="append",help="enter at least one output number",)
    parser.add_argument("-s","--step", type = int,default =1, help="enter at least one output number",)
    parser.add_argument("-p","--proj", type  = int, action="append", help="enter at least one output number",)        
    parser.add_argument("-ncols", type  = int,  help="enter at least one output number",)        
    parser.add_argument("-nrows", type  = int,  help="enter at least one output number",)        
    parser.add_argument("-midpoint", 	help = "Make plots up to output number n2. If None, sets n2 = n1.", type = int, default = None)
    parser.add_argument("--hist", action     = "store_true",default=False, help="enter at least one output number",)        
    parser.add_argument("-nolog", action     = "store_false",default=True, help="enter at least one output number",)        
    parser.add_argument("-average", action   = "store_false",default=True, help="enter at least one output number",)        
    parser.add_argument("-path",default="./", help="enter at least one output number",)        
    parser.add_argument("-multitime",default=None, help="list of snapshots",type =int, nargs="+")        
    parser.add_argument("-name",default=None, help="enter at least one output number",)        
    parser.add_argument("-impage",default=0, type = int,help="enter at least one output number",)        
    parser.add_argument("-rx", "--rotatex", action="store_true", default=False, help="Disable log scale for both axes")
    parser.add_argument("-ry", "--rotatey", action="store_true", default=False, help="Disable log scale for both axes")
    parser.add_argument("-rz", "--rotatez", action="store_true", default=False, help="Disable log scale for both axes")
    parser.add_argument("--timestamp", action="store_true", default=False, help="Disable log scale for both axes")
    args = parser.parse_args(subc_args, fix_n1n2=True)
    parser.set_outputdir(args.path+"/myplots/")

    dd = FastPlotter(args)
    dd.clim   = [
                 [2e-3 , 2e1 ], [2e-3, 2e1  ],
                 [2e2  , 1e8 ], [2e2 , 1e8  ],
                 [1e-9 , 1e-4], [1e-9, 1e-4 ],
                 [-2   , 5   ], [-2  , 5    ]
                 ]
    dd.run()