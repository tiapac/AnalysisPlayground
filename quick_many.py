import matplotlib.pyplot as plt
from multiprocessing import Pool
import os
import sys
import numpy as np
import utils.utils as utils
import yt
from yt.data_objects.level_sets.api import *
import analysis_playground.analysis_tools as an 
from matplotlib.colors import LogNorm

yt.funcs.matplotlib_style_context("dark_background")
import time as TIME
#initfield(MHD=False)####NEEDED TO INITIALIZE NEW FIELS
import argparse
from analysis_playground.constants import factG_in_cgs
plt.rcParams['axes.facecolor'] = 'black'
#
def add_pscalar(i):
	def _scalar(field, data):
		return (
			  data["gas", "density"]
			* data["ramses", "hydro_scalar_%02d"%i]
		)
	yt.add_field(
    name=("gas", "dscalar_%02d"%i),
    function=_scalar,
    sampling_type="local",
    #units="dyne/cm**2",
	)
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
 

parser = argparse.ArgumentParser(prog="Quickplot")
parser.add_argument("n1", 	help="enter at least one output number", type=int)
parser.add_argument("-n2", 	help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None)
#parser.add_argument("quantity",	help ="Quantity to analyse.", default = "density"   , type = str)
parser.add_argument("-nl","--nolog",  help ="plot linear variable."  , action  = "store_true", default=False)
parser.add_argument("-b","--box",	  help ="Box size to be plotted.", type = int)
parser.add_argument("-ax","--axis",	  help ="Set the axis to slice/project on.", default="x",type = str)
parser.add_argument("-p","--particles",help ="Add particles to the plot.", action  = "store_true", default=False)
parser.add_argument("-sb","--stream_B",help ="Add magnetic field streamlines.", action  = "store_true", default=False)
parser.add_argument("-sv","--stream_V",help ="Add velocity field streamlines.", action  = "store_true", default=False)
parser.add_argument("-pj","--projection",help ="Do a projection instead of a slice.", action  = "store_true", default=False)
parser.add_argument("--weight",help ="Weight to be used in the projection, if weighted.", default=None,type = str)
parser.add_argument("--path",help ="path to the directory containing outputs.", default="./",type = str)
parser.add_argument("--units",help ="Units of the colormap.", default=None,type = str)

parser.add_argument("--grid",help ="Plot the grid.", action  = "store_true", default=False)
parser.add_argument("-c","--center",help ="Center of the plot, three cooridinates between 0.0 <= x_i <= 1.0] are needed ", default =(0.5,0.5,0.5), nargs="+")
parser.add_argument("-d","--depth", help ="Depth of the projection along the --axis ", default = None,type=float)
parser.add_argument("--vlim",help ="Limits on the colormap.", default = None, type=float,      nargs=2)
parser.add_argument("--cmap",help ="Colormap to use.", default="jet",type = str)
parser.add_argument("--fields",help ="Show field and exit.", default=False,action  = "store_true")
parser.add_argument("-dir", help = "Name of the director where to store the profiles.", default = "./quikplots/data")
parser.add_argument("-tuo","--time_units_offset", help = "Time unit and offset in the same units.", default =["Myr",0.0] , nargs=2)
parser.add_argument("-xb","--xray_band", help = "minimum and maximum ernergy for the xray emission.", default =[0.5,3.0] , type= float, nargs=2)
parser.add_argument("--xray_source", help = "Model to compute xray spectrums.", default ="cloudy")#"spex" )
parser.add_argument("--arepo", help = "Analyse arepo files.", action="store_true", default=False )
parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)
parser.add_argument("-cpus", 	help="step in the data. Use only if n2 is used.", type=int, default=None)

parser.add_argument("-bcode","--box_code_units", 	help="use code units instade of image coordinates.", action  = "store_true", default=False)
parser.add_argument("-clumps", 	help="Run the clump finder.", action  = "store_true", default=False)
parser.add_argument("-clump_setting", 	help="Run the clump finder.", nargs=1)
parser.add_argument("-info", 	help="Print yt info.", action  = "store_true", default = False)
parser.add_argument("-scale", 	help="Annotate scale.", action  = "store_true", default=False)
parser.add_argument("-filter", 	help="Filter with a gaussian beam.", default=None, nargs=2, type=int)



#slc.annotate_arrow((0.5, 0.5, 0.5), length=0.06, color="blue")


args = parser.parse_args()
if args.cpus is not None:
    from multiprocessing import Pool
    from itertools import repeat
if not args.info: yt.set_log_level(50)### DO NOT SHOW THE INFO IN THE STDOUT

if not args.arepo:
	try:
		from mylibraries.addedfields import initfield
		initfield_exitsts = True
		initfield(MHD=True)
	except Exception as e :
		print(e, "File 'mylibraries.addedfields' not found.") 
		initfield_exitsts = False
args.time_units_offset[1]=float(args.time_units_offset[1])
args.center = np.array(args.center,dtype=float)
args.depth	= (args.depth,"pc")
if args.n2 is None: args.n2=int(args.n1)
if  args.n1 > args.n2: 
    n2tmp	= int(args.n2)
    args.n2 = int(args.n1)
    args.n1 = n2tmp
#quant_name = args.quantity
#weight_name= args.weight
#
 
#if "xray" in args.quantity: 
#    import pyxsim
#    source_model = pyxsim.CIESourceModel(args.xray_source,
#                                        emin   = 0.05, emax = 11., 
#                                        #kT_min = 0.025, kT_max = 64,
#                                        nbins  = 1000, 
#                                        #max_density=max_density,
#                                        Zmet   = 1.0,
#                                        binscale = "log" )
#    bandmin=args.xray_band[0]
#    bandmax=args.xray_band[1]
snaplist=range(args.n1, args.n2+1,args.step)                   
        
def look_data(i,args=args):
	outpath=os.getcwd() + "/"#"./"
	try:
		os.makedirs(args.dir)
	except:
		pass
	#if not args.arepo: initfield(MHD=False) ####NEEDED TO INITIALIZE NEW FIELS
	outname = outpath+"output_%05d"%i if not args.arepo else "snap_%03d.hdf5"%i
	ds = an.analyse(outname)
	#dataset=ds.ds
	t = ds.time.to("kyr")#.value
	print("Dataset %i loaded"%i)
	#if "xray" in args.quantity: 
	#	xray_fields  = source_model.make_source_fields(   dataset, bandmin,bandmax  )		
	#	print("Generated the X-ray fields:")
	#	print(xray_fields)
	#	#quant = ("gas", args.quantity+"_%s_%s_keV"%(str(bandmin),str(bandmax)))
	ds.init_sp()
	sp=ds.sp	
	vcut = 5.0e5
	tcut = 2e4
	cut_vel  = (sp["gas","velocity_magnitude"].in_cgs()).v > vcut
	field_tot = ["thermal_energy", "kinetic_energy", "magnetic_energy", "volume", "mass"]
	field_ave = ["density", "temperature", "pressure",
				"velocity_magnitude", "magnetic_field_magnitude", "magnetic_field_xz"]
	
	#regions = {"shell":cut_temp_shell, "bubble":cut_temp_bubble }
	print("Cutter data loaded.")
	def regions(name):
		if name=="bubble":
			result = cut_temp_bubble = (sp["gas","temperature"].in_cgs()).v > tcut
		else: 
			result = cut_temp_shell  = (sp["gas","temperature"].in_cgs()).v <= tcut
		return result 
	for name in ["shell", "bubble"]:
		weight = (sp["gas","mass"].to("Msun")).v
		cut_array = regions(name)
		cond = np.logical_and(cut_vel, cut_array)
		total_values = []		
		average_values = []
		std_values = []
		min_values = []
		max_values = []		
		for field in field_tot:
			#total
			tmp = (sp["gas",field].in_cgs()).v
			tmp = tmp[cond]
			total_values.append(tmp.sum())
			print("Done", field)
		weight = weight[cond]
		for field in field_ave:
			tmp = (sp["gas",field].in_cgs()).v
			tmp = tmp[cond]
			#average
			try:
				tmp_avg, tmp_std = utils.weighted_avg_and_std(tmp, weight)
			except: 
				tmp_avg=0.0
				tmp_std=0.0
			average_values.append(tmp_avg)
			std_values.append(tmp_std)
			#min and max
			try:
				min_values.append(tmp.min())
				max_values.append(tmp.max())
			except: 
				min_values.append(0.0)
				max_values.append(0.0)
			print("Done", field)
		del weight, cut_array,  cond, tmp
		total_values.append(t)
		average_values.append(t)
		std_values.append(t)
		min_values.append(t)
		max_values.append(t)
		total_values.append(i)
		average_values.append(i)
		std_values.append(i)
		min_values.append(i)
		max_values.append(i)
		total_values 	= np.array(total_values)
		average_values 	= np.array(average_values)
		std_values 		= np.array(std_values)
		min_values 		= np.array(min_values)
		max_values 		= np.array(max_values)
  
		with open(args.dir+"/%s_total.txt"%name, mode="+a") as FILE:
		    np.savetxt(FILE, total_values.reshape(1, total_values.shape[0]), fmt=' %.8e'*(len(total_values)-1)+" %d")
		with open(args.dir+"/%s_averages.txt"%name, mode="+a") as FILE:
		    np.savetxt(FILE, average_values.reshape(1, average_values.shape[0]), fmt=' %.8e'*(len(average_values)-1)+" %d")
		with open(args.dir+"/%s_std.txt"%name, mode="+a") as FILE:
		    np.savetxt(FILE, std_values.reshape(1, std_values.shape[0]), fmt=' %.8e'*(len(average_values)-1)+" %d")
		with open(args.dir+"/%s_min.txt"%name, mode="+a") as FILE:
		    np.savetxt(FILE, min_values.reshape(1, min_values.shape[0]), fmt=' %.8e'*(len(average_values)-1)+" %d")
		with open(args.dir+"/%s_max.txt"%name, mode="+a") as FILE:
		    np.savetxt(FILE, max_values.reshape(1, max_values.shape[0]), fmt=' %.8e'*(len(average_values)-1)+" %d")
			
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
