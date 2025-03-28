from utils.myparser import myparser
from yt.data_objects.level_sets.api import *

def MakeHist(parser = myparser("Histograms"), subc_args = None ):
	import matplotlib.pyplot as plt
	from multiprocessing import Pool
	import os
	import sys
	import numpy as np
	import utils.utils as utils
	import yt
	import analysis_playground.analysis_tools as an 
	from matplotlib.colors import LogNorm

	# yt.funcs.matplotlib_style_context("dark_background")
	import time as TIME
	#initfield(MHD=False)####NEEDED TO INITIALIZE NEW FIELS
	import argparse
	from analysis_playground.constants import factG_in_cgs
	#plt.rcParams['axes.facecolor'] = 'black'

	

	# parser = argparse.ArgumentParser(prog="Quick_hist")
	parser.add_argument("n1", 	help="enter at least one output number", type=int)
	parser.add_argument("quantity", 	help="enter at least one output number", type=str)
	parser.add_argument("-n2", 	help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None)
	parser.add_argument("-bins",  help = "Number of bins to compute averages"   , type=int, default = 256)
	parser.add_argument("-d","--density",  help ="plot linear variable on x."  , action  = "store_true", default=False)

	parser.add_argument("-nlx","--nologx",  help ="plot linear variable on x."  , action  = "store_true", default=False)
	parser.add_argument("-nly","--nology",  help ="plot linear variable on y."  , action  = "store_true", default=False)
	parser.add_argument("--weight",    help ="Weight to be used in the projection, if weighted.", default=None, type = str)
	parser.add_argument("-c","--color",help ="Color of the histogram.",    default = "midnightblue",    nargs="+")
	parser.add_argument("-type",help ="Hist type",choices=['bar', 'barstacked', 'step', 'stepfilled'],    default = "bar")
	parser.add_argument("-ec",help ="Hist edge color", default = None)
	parser.add_argument("-dir", help = "Name of the director where to store the profiles.", default = "./quikplots")
	parser.add_argument("-cpus", 	help="step in the data. Use only if n2 is used.", type=int, default=None)

	parser.add_argument("-paths",      help ="Paths to the dataset/s.", default = "./",       nargs=1,)
	parser.add_argument("-yl","--ylabel",help ="Label/s for the hist/s.", default = None,       nargs=1,)
	parser.add_argument("-xl","--xlabel",help ="Label/s for the hist/s.", default = None,       nargs=1,)

	#parser.add_argument("-time","--time_label",help ="Use time as label",  action  = "store_true", default = False)

	parser.add_argument("-uw",help ="units on the weight", default = 1.0,   type=float,    nargs=1)
	parser.add_argument("-ux",help ="units on the x axis", default = 1.0,   type=float,    nargs=1)

	parser.add_argument("--ylim",help ="Limits on the y axis", default = None,    type=float,   nargs=2)
	parser.add_argument("--xlim",help ="Limits on the x axis", default = None,    type=float,   nargs=2)
	parser.add_argument("-abs","--absolute",help ="Use absolute values for the quantity.",  action  = "store_true", default = False)
	parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)
	


	args = parser.parse_args(subc_args)
	if args.cpus is not None:
		from multiprocessing import Pool
		from itertools import repeat
	#if not args.info: yt.set_log_level(50)### DO NOT SHOW THE INFO IN THE STDOUT

	if args.n2 is None: args.n2=int(args.n1)
	if  args.n1 > args.n2: 
		n2tmp	= int(args.n2)
		args.n2 = int(args.n1)
		args.n1 = n2tmp

	snaplist=range(args.n1, args.n2,args.step)                   
	# def check_dir(): 
	args.dir += "/" +"quick_histograms"
	args.dir += "/" +args.quantity+"/"
	try:
		os.makedirs(args.dir)
	except OSError as error:
		print(error)  
	
		
		
		
		
	def look_data(i,args=args):
		outpath=os.getcwd() + "/"#"./"
		try:
			os.makedirs(args.dir)
		except:
			pass
		#if not args.arepo: initfield(MHD=False) ####NEEDED TO INITIALIZE NEW FIELS
		outname = outpath+"output_%05d"%i 
		ds = yt.load(outname)
		t  = ds.current_time.to("kyr")#.value
		print("Dataset %i loaded"%i)
		sp = ds.all_data()
		quant = sp["gas", args.quantity].v/args.ux
		if args.absolute: quant = abs(quant)
		if args.weight is not None: 
			weight = sp["gas", args.weight].v/args.uw
		else: 
			weight = None
		if not args.nologx: quant = np.log10(quant) 
		
		fig = plt.figure(1)
		ax = fig.add_subplot(111)
		if weight is None: 
			ax.hist(quant, bins=args.bins, ec=args.ec, log=not args.nology, color=args.color, density=args.density )
		else:
			ax.hist(quant, bins=args.bins, weights=weight, ec=args.ec, log=not args.nology, color=args.color, density=args.density )
		if args.xlabel is not None: ax.set_xlabel(r"%s"%args.xlabel)
		if args.ylabel is not None: ax.set_ylabel(r"%s"%args.ylabel)
		if args.xlim is not None: ax.set_xlim(args.xlim[0],args.xlim[1])
		if args.ylim is not None: ax.set_ylim(args.ylim[0],args.ylim[1])
		fig.savefig(args.dir+"/hist_%i_%s"%(i,args.quantity))
	
	
		return 
		

		



	if args.cpus is not None:
		with Pool(processes=args.cpus)as p:
			results=p.starmap_async(look_data, zip(snaplist, repeat(args)) )
			results.get()###NEEDED TO OBTAIN THE RESULTS IN ASYNC
			p.close()
	else:
		for i in snaplist:
			look_data(i,args)



	#s.annotate_cell_edges()
	print("Done.")
