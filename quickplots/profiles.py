import analysis_playground.analysis_tools as an
import matplotlib.pyplot as plt
import os 
import numpy as np
import argparse
 
parser = argparse.ArgumentParser(prog="Quickprofile")
parser.add_argument("n1", 	help="enter at least one output number", type=int, nargs="+")
parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)
parser.add_argument("-n2",  help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None,nargs="+")
parser.add_argument("-dir", help = "Name of the director where to store the profiles.", default = "./quikplots")
parser.add_argument("-o","--output", 	help="Name of the ouput file in fig.savefig(...).", default = None)
parser.add_argument("--path",help ="path to the directory containing outputs.", default=None,type = str)

parser.add_argument("x",help ="Quantity on x axis."   , type = str)
parser.add_argument("y",help ="Quantity on y axis."   , type = str)

parser.add_argument("-bins",  help = "Number of bins to compute averages"   , type=int, default = 128)
parser.add_argument("-nlx","--nologx",  help ="plot linear variable on x."  , action  = "store_true", default=False)
parser.add_argument("-nly","--nology",  help ="plot linear variable on y."  , action  = "store_true", default=False)

parser.add_argument("--weight",    help ="Weight to be used in the projection, if weighted.", default="mass", type = str)
parser.add_argument("-c","--color",help ="Color of the line/s.",    default = ["black"],    nargs="+")
parser.add_argument("-cmap",help ="Map the color on a color map. TO be only used for multiple ouputs, not different simulations.",    default = None)

parser.add_argument("-paths",      help ="Paths to the dataset/s.", default = ["./"],       nargs="+",)
parser.add_argument("-l","--label",help ="Label/s for the line/s.", default = None,       nargs="+",)
parser.add_argument("-time","--time_label",help ="Use time as label",  action  = "store_true", default = False)

parser.add_argument("--linestyle", help ="stile of the line/s.",  default = ["-"],        nargs="+",)
parser.add_argument("--single_plot",help ="Decide if the plots should be drwan all in the same figure.",  action  = "store_true", default = False)
parser.add_argument("--additional_fields", 	help="Additional fields to plot.", type=str, nargs="+")


parser.add_argument("--ylim",help ="Limits on the y axis", default = None,       nargs=2)
parser.add_argument("--xlim",help ="Limits on the x axis", default = None,       nargs=2)


args = parser.parse_args()
if args.additional_fields is not None and len(args.linestyle)<len(args.additional_fields):
    args.linestyle=args.linestyle*(len(args.additional_fields)+1)

nfiles = len(args.n1)
if args.xlim is not None: args.xlim = np.array(args.xlim,dtype=float)
if args.ylim is not None: args.ylim = np.array(args.ylim,dtype=float)
def check_dir(): 
    try:
        os.mkdir(args.dir)
    except OSError as error:
        print(error) 
    args.dir += "/" +"profiles"
    try:
        os.mkdir(args.dir)
    except OSError as error:
        print(error) 
    try:
        os.mkdir(args.dir)
    except OSError as error:
        print(error) 
        
    args.dir += "/" +args.x

    try:
        os.mkdir(args.dir)
    except OSError as error:
        print(error) 
        
    args.dir += "/" +args.y

    try:
        os.mkdir(args.dir)
    except OSError as error:
        print(error) 
if args.n2 is None: args.n2 = args.n1


if args.cmap is not None: 
    cmap=plt.get_cmap(args.cmap)
    args.color = cmap(np.linspace(0,1,
                                    (max(args.n2)-min(args.n1))//(args.step)+1))

if args.single_plot:
            fig = plt.figure(0)
            ax = fig.add_subplot(111)
            
first=True         
for j in range(0,nfiles):
    k = 0

    for i in range(args.n1[j],args.n2[j]+1,args.step):
        #print(i,k, len(args.color),len(range(args.n1[j],args.n2[j]+1,args.step)))
        if not args.single_plot:
            fig = plt.figure(i)
            ax = fig.add_subplot(111)
        ds = an.analyse(outpath=args.paths[j] ,outnumb = i)
        if args.label is None:
            label = None
            if args.time_label: label = " t = %.3E Myr"%ds.time.to("Myr").v
              #else 
        else: 
            label = args.label[j]
            if args.time_label: label += " t = %.3E Myr"%ds.time.to("Myr").v
        print("Plotting",i)
        #print(args.linestyle)
        if args.additional_fields is not None:
            for l,yy in enumerate(args.additional_fields):
                ax,fig = ds.oneplot(args.x,yy,
                        ax = ax, fig = fig,
                        binned       = True, 
                        w            = args.weight,
                        bins         = args.bins,
                        color        = args.color[j if args.cmap is None else k],
                        logx         = not args.nologx,
                        logy         = not args.nology,
                        #label        = label,
                        linestyle    = args.linestyle[l+1]
                        )
        ax,fig = ds.oneplot(args.x,args.y,
                   ax = ax, fig = fig,
                   binned       = True, 
                   w            = args.weight,
                   bins         = args.bins,
                   color        = args.color[j] if args.cmap is None else args.color[k],
                   logx         = not args.nologx,
                   logy         = not args.nology,
                   label        = label,
                   linestyle    = args.linestyle[j]
                   )
        
        

        if args.ylim is not None: ax.set_ylim(args.ylim[0], args.ylim[1])
        if args.xlim is not None: ax.set_xlim(args.xlim[0], args.xlim[1])
        if first: check_dir()
        name = args.dir + "/" + args.x+"_"+args.y + "_%05d"%i
        if args.output is not None: name = args.output
        if not args.single_plot: 
            fig.savefig(name)
            plt.close(fig)
        else:
            try:
                os.mkdir(args.dir + "/tmp")
            except OSError as error:
                print(error)    
            name = args.dir + "/tmp/" + args.x+"_"+args.y + "_%05d"%i
            fig.savefig(name)
             
        k += 1
        first=False
name = args.dir +"/" + args.x+"_"+args.y  +"_total_%05d-%05d"%(min(args.n1),max(args.n2))
if args.output is not None: name = args.output
if args.single_plot: fig.savefig(name)
        