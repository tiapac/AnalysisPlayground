from utils.myparser import myparser

import analysis_playground.analysis_tools as an
import matplotlib.pyplot as plt
import os
import numpy as np
import paper_settings as pg
import constants as cc

def make_a_phase(self, i):
        #print("Doing {i}")
        print(f"Plotting phase of snapshot {i}")
        if self.args.weight is not None:
            weigthtype = "gas" if "hydro" not in self.args.weight else "ramses"
        else:
            weigthtype = None

        fig = plt.figure(0)
        ax = fig.add_subplot(111)
        fig.set_size_inches(self.paper.onecol, self.paper.onecol)
        ds = an.analyse(outnumb=i)
        print(f"Plotting phase of snapshot {i} at time {ds.time.to('kyr').value}")
        ax, fig = ds.phase(
            self.args.x, self.args.y,
            ax=ax, fig=fig,
            w=self.args.weight,
            weigthtype=weigthtype,
            bins=self.args.bins,
            logx=not self.args.nologx,
            logy=not self.args.nology,
            logw=not self.args.nologw,
            sumist=self.args.sum,
            shading='gouraud',
            only_scalar=self.args.only_scalar,
            absx=self.args.absolute_x,
            absy=self.args.absolute_y,
            cmap=self.args.cmap,
            **self.phaseargsD
        )
        if self.args.ylim is not None:
            ax.set_ylim(self.args.ylim[0], self.args.ylim[1])
        if self.args.xlim is not None:
            ax.set_xlim(self.args.xlim[0], self.args.xlim[1])

        if self.args.lines is not None:
            for lines in self.args.lines:
                a, b, xxmin, xxmax, color = map(float, lines[:4]) + [lines[4]]
                xarray = np.linspace(xxmin, xxmax, num=1000)
                ax.plot(xarray, a + (xarray - xxmin) * b, color=color)

        name = f"{self.args.dir}/phase_{self.args.x}_{self.args.y}_{self.wname}_{i:05d}.png"
        if self.args.output is not None:
            name = self.args.output

        dpi = self.paper.figure.dpi
        if self.args.curve is not None:
            c_path, x_idx, y_idx = self.args.curve
            curve_data = np.loadtxt(c_path)
            cc_x, cc_y = curve_data[:, int(x_idx)], curve_data[:, int(y_idx)]
            ax.plot(cc_x, cc_y, color="black", lw=self.paper.figure.linewidth, alpha=0.6)
            gamma = 5. / 3.
            n1 = 0.7
            T1 = 6000.0
            rho1 = n1 * cc.mH * cc.mu_mol
            P1 = T1 * n1 * cc.kB
            e1 = P1 / ((gamma - 1) * rho1)
            rho2 = np.linspace(rho1, 5 * rho1, 500)
            P2 = P1 * ((gamma + 1) / (gamma - 1)) * (rho2 / rho1) - (2 * P1 / (gamma - 1))
            ax.plot(rho2, P2)

        fig.set_size_inches(*self.paper.psize[self.args.impage])
        fig.tight_layout()
        fig.savefig(name, dpi=dpi)
        plt.close(fig)
        print("Done")
class PhasePlot:
    def __init__(self, args):
        self.args = args
        self.paper = pg.PaperSettings()
        self.phaseargsD = self.parse_phaseargs(args.phaseargs)
        self.snaplist = range(args.n1, args.n2 + 1, args.step)
        self.start = args.n1
        self.end = args.n2
        self.step = args.step
        self.wname = args.weight if args.weight is not None else "num"
        self.setup_directories()
        plt.rcParams.update({"text.usetex": True, **self.paper.mplfigfontstyle})
        

    def parse_phaseargs(self, phaseargs):
        phaseargsD = {}
        for idx, arg in enumerate(phaseargs[:-1]):
            if arg.startswith("--") or arg.startswith("-"):
                key = arg
                if not phaseargs[idx + 1].startswith("-") or phaseargs[idx + 1].startswith("--"):
                    phaseargsD[key.replace("-", "")] = phaseargs[idx + 1]
                else:
                    phaseargsD[key.replace("-", "")] = [phaseargs[idx + 1]]
                    tempidx = idx + 2
                    while not (phaseargsD[key].startswith("-") or phaseargsD[key].startswith("--")) or tempidx + 1 == len(phaseargs):
                        phaseargsD[key] += phaseargs[idx]
                        tempidx += 1
        return phaseargsD

    def setup_directories(self):
        self.args.dir += "/" + "phase_plots"
        self.args.dir += "/" + self.args.x
        self.args.dir += "/" + self.args.y
        self.args.dir += "/" + self.wname
        try:
            os.makedirs(self.args.dir)
        except:
            pass

    

    def run(self):
        snaplist = range(self.start, self.end+1, self.step)
        if self.args.cpus is not None:
            from itertools import repeat
            from multiprocessing import Pool
            
            with Pool(processes=self.args.cpus) as p:
                results = p.starmap_async(make_a_phase, zip(repeat(self),snaplist))
                results.get()
                p.close()
        else:
            for i in snaplist:
                
                make_a_phase(self,i)

def MakePhase(parser = myparser(prog = "Fastplosts"), subc_args = None):
    parser.add_argument("n1", help="enter at least one output number", type=int)
    parser.add_argument("x", help="Quantity on x axis.", type=str)
    parser.add_argument("y", help="Quantity on y axis.", type=str)
    parser.add_argument("-s", "--step", help="step in the data. Use only if n2 is used.", type=int, default=1)
    parser.add_argument("-n2", help="Make plots up to output number n2. If None, sets n2 = n1.", type=int, default=None)
    parser.add_argument("-dir", help="Name of the directory where to store the profiles.", default="./quikplots")
    parser.add_argument("--output", help="Name of the output file in fig.savefig(...).", default=None)
    parser.add_argument("-bins", help="Number of bins to compute averages for the two quantity. Can be just one.", type=int, default=[256], nargs="+")
    parser.add_argument("-nlx", "--nologx", help="plot linear variable on x.", action="store_true", default=False)
    parser.add_argument("-nly", "--nology", help="plot linear variable on y.", action="store_true", default=False)
    parser.add_argument("-nlw", "--nologw", help="plot linear variable on the weight.", action="store_true", default=False)
    parser.add_argument("--sum", help="Plot the sum rather than the averages in each bin.", action="store_true", default=False)
    parser.add_argument("-w", "--weight", help="Weight to be used in the projection, if weighted.", default=None)
    parser.add_argument("-cmap", help="Colormap to be used.", default="cividis")
    parser.add_argument("-cpus", help="step in the data. Use only if n2 is used.", type=int, default=None)
    parser.add_argument("--ylim", help="Limits on the y axis", default=None, nargs=2)
    parser.add_argument("--xlim", help="Limits on the x axis", default=None, nargs=2)
    parser.add_argument("--only_scalar", help="Only consider cells with passive scalar 2 concentration > 1e-15", action="store_true", default=False)
    parser.add_argument("-absx", "--absolute_x", help="Use absolute values for the quantity.", action="store_true", default=False)
    parser.add_argument("-absy", "--absolute_y", help="Use absolute values for the quantity.", action="store_true", default=False)
    parser.add_argument("--lines", help="a and b, x-min and -max, color to draw n straight line", default=None, nargs=5, action="append")
    parser.add_argument("-curve", help="plot curve from data.", nargs=3, default=None)
    parser.add_argument("-impage", help="plot curve from data.", default=0, type=int)

    args, phaseargs = parser.parse_known_args(subc_args)
    args.phaseargs = phaseargs
    
    if args.n2 is None:
        args.n2 = args.n1
    if args.n1 > args.n2:
        args.n1, args.n2 = args.n2, args.n1
    if len(args.bins) == 1:
        args.bins.append(args.bins[0])
    if args.xlim is not None:
        args.xlim = np.array(args.xlim, dtype=float)
    if args.ylim is not None:
        args.ylim = np.array(args.ylim, dtype=float)
    
    phase_plot = PhasePlot(args)
    phase_plot.run()

