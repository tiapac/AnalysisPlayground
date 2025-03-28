import matplotlib.pyplot as plt
from multiprocessing import Pool
import os
import sys
import numpy as np
import utils.utils as utils
import yt

import time as TIME
#initfield(MHD=False)####NEEDED TO INITIALIZE NEW FIELS
import argparse
yt.set_log_level(50)### DO NOT SHOW THE INFO IN THE STDOUT
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

parser = argparse.ArgumentParser(prog="Quickplot")
parser.add_argument("n1", 	help="enter at least one output number", type=int)
parser.add_argument("quantity",	help ="Quantity to analyse.", default = "density"   , type = str)
parser.add_argument("-nl","--nolog",  help ="plot linear variable."  , action  = "store_true", default=False)
parser.add_argument("-b","--box",	  help ="Box size to be plotted.", type = int)
parser.add_argument("-ax","--axis",	  help ="Set the axis to slice/project on.", default="x",type = str)
parser.add_argument("-p","--particles",help ="Add particles to the plot.", action  = "store_true", default=False)
parser.add_argument("-sb","--stream_B",help ="Add magnetic field streamlines.", action  = "store_true", default=False)
parser.add_argument("-sv","--stream_V",help ="Add velocity field streamlines.", action  = "store_true", default=False)
#parser.add_argument("-pj","--projection",help ="Do a projection instead of a slice.", action  = "store_true", default=False)
parser.add_argument("--weight",help ="Weight to be used in the projection, if weighted.", default=None,type = str)
parser.add_argument("--grid",help ="Plot the grid.", action  = "store_true", default=False)
parser.add_argument("-c","--center",help ="Center of the plot, three cooridinates between 0.0 <= x_i <= 1.0] are needed ", default =(0.5,0.5,0.5), nargs="+")
parser.add_argument("-d","--depth", help ="Depth of the projection along the --axis, start and finish coordinates.", default = [0.,1.] ,type=float, nargs=2)
parser.add_argument("--vlim",help ="Limits on the colormap.", default = None, type=float,      nargs=2)
parser.add_argument("--cmap",help ="Colormap to use.", default="jet",type = str)
parser.add_argument("--fields",help ="Show field and exit.", default=False,action  = "store_true")
parser.add_argument("-dir", help = "Name of the director where to store the profiles.", default = "./quikplots/maps")
parser.add_argument("-tuo","--time_units_offset", help = "Time unit and offset in the same units.", default =["Myr",0.0] , nargs=2)
parser.add_argument("-xb","--xray_band", help = "minimum and maximum ernergy for the xray emission.", default =[0.5,3.0] , type= float, nargs=2)
parser.add_argument("--xray_source", help = "Model to compute xray spectrums.", default ="spex" )
parser.add_argument("--arepo", help = "Analyse arepo files.", action="store_true", default=False )
parser.add_argument("-s","--step", 	help="steps to cover the --depth.", type=int, default=1)
parser.add_argument("-bcode","--box_code_units", 	help="use code units instead of image coordinates.", action  = "store_true", default=False)

args = parser.parse_args()
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
#args.depth	= (args.depth,"pc")
quant_name = args.quantity
weight_name= args.weight
quant= ('gas', quant_name) if "hydro" not in quant_name  else ("ramses", quant_name)
args.step=args.step-1
if args.weight is not None: 
    weight_field = ('gas', weight_name) if "hydro" not in weight_name  else ("ramses",weight_name)
else: 
    weight_field = None
if "dscalar"in quant_name:  
    for i in range(1,2+1): 
        add_pscalar(i)
 
if "xray" in args.quantity: 
    import pyxsim
    source_model = pyxsim.CIESourceModel(args.xray_source,
                                        emin   = 0.05, emax = 11., 
                                        #kT_min = 0.025, kT_max = 64,
                                        nbins  = 1000, 
                                        #max_density=max_density,
                                        Zmet   = 1.0,
                                        binscale = "log" )
    bandmin=args.xray_band[0]
    bandmax=args.xray_band[1]
                   
    
told=0.
outpath=os.getcwd() + "/"#"./"
outname = outpath+"output_%05d"%args.n1 if not args.arepo else "snap_%03d.hdf5"%args.n1
ds = yt.load(outname,default_species_fields='ionized')
if "xray" in args.quantity: 
		xray_fields  = source_model.make_source_fields(   ds, bandmin,bandmax  )
		#xray_fields += source_model.make_intensity_fields(ds, bandmin,bandmax, dist=(0.0,"pc")  )
		
		print("Generated the X-ray fields:")
		print(xray_fields)
		quant = ("gas", args.quantity+"_%s_%s_keV"%(str(bandmin),str(bandmax)))

if args.fields: 
	for f in ds.derived_field_list:
		print(f)
	quit()
boxlen = float(ds.domain_width[0].to("pc").value)
if args.box is None:
	width = float(ds.domain_width[0].to("pc").value)
else: 
	width = args.box
	boxlen = float(ds.domain_width[0].to("pc").value)

if args.box_code_units:   
	dx = ( (args.depth[1]-args.depth[0])/ boxlen ) / args.step 
 
else:
    print(args.depth[1],args.depth[0])
    dx = (args.depth[1]-args.depth[0]) / args.step 


if args.box_code_units:
		center = abs( np.array(args.center) + boxlen * 0.5 ) / boxlen
else: 
		center= (args.center).copy()
    		
ax_idx = 0
for ax in ["x","y","z"]:
    if args.axis == ax: break
    ax_idx += 1
    
for i in range( 0 , args.step + 1 ):
	
	t = ds.current_time.to( "Myr" ) 
	center[ax_idx] = i * dx if i>0 else args.depth[0]  # (args.center).copy()
	if i * dx ==1.:
		center[ax_idx] =  0.999  # (args.center).copy()
	
	if args.box_code_units:
		center = abs( np.array(args.center) + boxlen * 0.5 ) / boxlen
	
  
	#if not args.projection:
	print( boxlen, center)
	s  = yt.SlicePlot(ds, args.axis, quant, width=(width,"pc"), center=center) #, center=center) #, no_ghost=False)### do the slice plto
  
	if "velocity" in quant_name and not "div" in quant_name  :s.set_unit(quant, "km/s")
 
	s.set_cmap(field="all", cmap=args.cmap)

	if  args.nolog:		
		s.set_log(field = quant, log = False)
  
	if args.grid: s.annotate_grids()
 
	s.annotate_timestamp(time_offset=(args.time_units_offset[1],args.time_units_offset[0]), time_unit=args.time_units_offset[0])
	s.annotate_text((0.8, 0.8),
                 r"$%s_0=%.2f \rm pc $"%(args.axis, center[ax_idx]*boxlen  ),
                 coord_system="axis")
	axses = utils.find_axes(args.axis)

	if args.axis == "y":
			s.swap_axes()
			axses.reverse()

	if args.stream_B:	
		s.annotate_streamlines(("gas","magnetic_field_"+axses[0]),("gas","magnetic_field_"+axses[1]),
                                   density=5,
                                   color="midnightblue",
                                   arrowsize=0.5,
                                   )
	if args.stream_V:
		s.annotate_streamlines(("gas","velocity_"+axses[0]),("gas","velocity_"+axses[1]),
	                               density=5,
	                               color="crimson",
	                               arrowsize=0.5,
	                               )
  
	if args.particles: s.annotate_particles(width = width, marker = "*", p_size = 100,)

	if args.vlim is not None: s.set_zlim(quant, args.vlim[0],args.vlim[1])

	plottype="mri_slices"  
	with_depth="mri_%05d"%args.n1
	s.save(args.dir+"/"+plottype+"/"+with_depth+"/"+quant_name+"/"+args.axis+"/%05d" %i)
	print(i, "Time=%.4f"%t)#,"Tmax=%.4f"%temp.min())


#s.annotate_cell_edges()
print("Done.")
