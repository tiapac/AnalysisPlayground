import analysis_playground.analysis_tools as an
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse
import yt
def _jeansLength(field, data):
    G = yt.YTArray(6.67430e-8, "cm**3/g/s**2") #G = factG_in_cgs*(yt.units.dyne*yt.units.cm/yt.units.s**2)
    cs =  data["gas", "sound_speed"].in_cgs()
    rho=  data["gas", "density"].in_cgs()
    return np.sqrt(np.pi * cs**2 / (G * rho))
    
yt.add_field(
    name=("gas", "jeans_length"),
    function=_jeansLength,
    sampling_type="local",
    #units="cm",
	)
 
parser = argparse.ArgumentParser(prog="Quickprofile")
parser.add_argument("n1", 	help="enter at least one output number", type=int, nargs="+")
parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)
parser.add_argument("-n2",  help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None,nargs="+")
parser.add_argument("-dir", help = "Name of the director where to store the profiles.", default = "./quikplots")
parser.add_argument("-o","--output", 	help="Name of the ouput file in fig.savefig(...).", default = None)

parser.add_argument("x",help ="Quantity on x axis."   , type = str)

parser.add_argument("--weight",    help ="Weight to be used in the projection, if weighted.", default="mass", type = str)
parser.add_argument("-paths",      help ="Paths to the dataset/s.", default = ["./"],       nargs="+",)

args = parser.parse_args()
if args.n2 is None: args.n2 = args.n1



nfiles = len(args.n1)
     
for j in range(0,nfiles):
    k = 0
    for i in range(args.n1[j],args.n2[j]+1,args.step):
        #print(i,k, len(args.color),len(range(args.n1[j],args.n2[j]+1,args.step)))
     
        ds = an.analyse(outpath=args.paths[j] ,outnumb = i)
        print(ds.q_mean(args.x))
