
import analysis_playground.analysis_tools as an
import numpy as np 
import constants as costs
import utils_3D
import vtk

import pyvista as pv
from pyvista import opacity_transfer_function as transfer 
import yt
#import transfer_functions as tff 
#import matplotlib.pyplot as plt
import sys
from utils.utils import flatten
import os
import argparse

#filename  = "./output_00025/"


parser = argparse.ArgumentParser(prog="Quickplot3D")
parser.add_argument("n1", 	help="enter at least one output number", type=int)
parser.add_argument("-n2", 	help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None)
parser.add_argument("-quantity",	help ="Quantity to analyse.", default = "density"   , type = str)
parser.add_argument("-nl","--nolog",  help ="plot linear variable."  , action  = "store_true", default=False)
parser.add_argument("-b","--box",	  help ="Box size to be plotted.", type = int, default=None)
#parser.add_argument("-ax","--axis",	  help ="Set the axis to slice/project on.", default="x",type = str)
parser.add_argument("-p","--particles",help ="Add particles to the plot.", action  = "store_true", default = False)
#parser.add_argument("-c","--center",help ="Center of the plot, three cooridinates between 0.0 <= x_i <= 1.0] are needed ", default =(0.5,0.5,0.5), nargs="+")
parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)

parser.add_argument("--cmap",help ="Colormap to use.", default="jet",type = str)
parser.add_argument("-dir", help = "Name of the director where to store the image.", default = "./quikplots/3D")
#parser.add_argument("-tuo","--time_units_offset", help = "Time unit and offset in the same units.", default =["Myr",0.0] , nargs=2)
#parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)
#parser.add_argument("-bcode","--box_code_units", 	help="use code units instade of image coordinates.", action  = "store_true", default=False)
#parser.add_argument("-clumps", 	help="Run the clump finder.", action  = "store_true", default=False)
#parser.add_argument("-clump_setting", 	help="Run the clump finder.", nargs=1)
parser.add_argument("-info", 	help="Print yt info.", action  = "store_true", default = False)
parser.add_argument("-sphere", 	help="Points as spheres.", action  = "store_true", default = False)
parser.add_argument("-gaussian", 	help="Points as gaussian spheres.", action  = "store_true", default = False)
parser.add_argument("-show", 	help="show plot.", action  = "store_true", default = False)

parser.add_argument("-load_data", 	help="Load the data.", action  = "store_true", default = True )
parser.add_argument("-load_grid", 	help="Load the grid.", action  = "store_true", default = True )
parser.add_argument("-save_grid", 	help="Save the grid.", action  = "store_true", default = True )
parser.add_argument("-save_data", 	help="Save the data.", action  = "store_true", default = False )
parser.add_argument("-screenshot", 	help="Save a screenshot. Set all off screen.", action  = "store_true", default = True )

parser.add_argument("-trhf","--treshold_field", 	help="Treshold field.", action  = "append", default = ["hydro_scalar_00"] )
parser.add_argument("-trhf_v","--treshold_field_value", help="--treshold field value.", action  = "append", default = [1e-15], nargs = 1 , type = float)
parser.add_argument("-trhf_o","--treshold_field_operation", help="--treshold field operation.", action  = "append", default = [">"] )
parser.add_argument("-diffuse", help="Diffuse light from points", default = 0.1, type = float )
parser.add_argument("-ps","--point_size", help="size of the points", default = 0.5, type = float )
parser.add_argument("-rescale", help="rescale the dataset n-times.", default = -1, type = int )
parser.add_argument("-cpus", 	help="step in the data. Use only if n2 is used.", type=int, default=None)

args = parser.parse_args()
dirs = (args.dir[:-1]) if (args.dir).endswith("/") else args.dir
dirs = dirs.split("/")
for folder in dirs:
    try:
        os.mkdir(folder)
    except OSError as error:
            print(error) 
if args.info: yt.set_log_level(0)### DO NOT SHOW THE INFO IN THE STDOUT
if args.cpus is not None:
    from multiprocessing import Pool
    from itertools import repeat

pv.global_theme.allow_empty_mesh = True
path                     = args.dir
quant                    = args.quantity#"xray_intensity_0.05_4.0_keV"#key+"_magnitude"
boxlen                   = args.box if args.box is None else [[-args.box/2.,args.box/2.]]*3
bounds                   = boxlen
axes                     = ["x","y","z"]#,["magnetic_field_magnitude","density","temperature"]
clim                     = None         # [1.e3,1.e4]#None#[max(qthreshold, cmin), cmax] #None
render_points_as_spheres = args.sphere

if not os.path.exists(path): os.mkdir(path)


#key  = "velocity"
load_data = True
load_grid = True
save_grid = True
save_data = False
#plotter = Plotter(lighting='none', off_screen=False)
istart    = args.n1 #int(sys.argv[1])
plotter   = pv.Plotter( off_screen = args.screenshot)
#plotter.camera.azimuth=0
NR_FRAMES   =  0        # len(ts)-istart
dtheta    =   0. if NR_FRAMES<=0 else  360./NR_FRAMES
poss      =   np.array((plotter.generate_orbital_path(factor=800.0, n_points=NR_FRAMES+1, viewup=(0,0,1), shift=800.0)).points, dtype=float)
plotter.camera.position=poss[0]
first   =   True
    
def tff( model = None, i = 32, ncolors=64, alpha=0):
        n    = 3
        t    = np.linspace(0, np.pi * n, ncolors)+np.pi/6
        func =  np.sin(3*t)*(1 - t / t.max())**( alpha)
        func = t * 0.5#np.sin(3*t)*(1 - t / t.max())**( alpha)
        func[func<0.] = 0.1
        
        opacity   = transfer(func, 256).astype(float)

        return  opacity
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
)
pv.OFF_SCREEN = True
if args.n2 is None: args.n2=int(args.n1)
def make_a_plot(i, args, camera_stuff=None):
    
    print("working on %i"%i)
    ds  = an.analyse(outnumb=i)
    time = ds.time.to("Myr")
    plotter=pv.Plotter(off_screen=not args.show)
    pltt = ds.plot3D(quant, render_points_as_spheres = True,
        axes            = axes,
        load_data       = load_data,
        save_data       = save_data,
        threshold_field = args.treshold_field,
        threshold       = args.treshold_field_value,
        operation       = args.treshold_field_operation,
        gaussian        = args.gaussian, 
        #log_axes=[True,True,True],
        #scale_axes=[1,1,1],
        clim            = clim,
        tf              = tff(alpha=0.5),#"linear" ,#tff(alpha=0.)*0.9, 
        #log=True ,
        point_size      = args.point_size,#.5,#0.8,#"adaptive",#1.5, 
        diffuse         = args.diffuse,   # 1.0,#0.1,# "adaptive",#0.0075,
        particles       = args.particles, 
        trajectories_from_file = False,
        trajectories_path="/fast/pacicco/ramses_inkTurb/ramses_ink/mypatch_mhdMi_star/mine/trajectories/", 
        density_phase   = axes != ["x","y","z"],
        plotter         = plotter, 
        centered        = True,
        rescale=args.rescale,
        )



    plotter.set_background("black")
    #if first:
    #    pltt.camera.is_set = True
    #    first              = False
    plotter.show_bounds(
                            bounds=bounds if bounds is None else flatten(bounds),
                            use_3d_text= False,
                            color      = "white",
                            bold       = False,
                            location   = "origin",#'outer',                       
                            ticks      = 'both',                       
                            n_xlabels  = 4,                        
                            n_ylabels  = 4,                        
                            n_zlabels  = 4,                        
                            xtitle     = "x",                       
                            ytitle     = "y",                      
                            ztitle     = "z",    
                            font_size  = 20,        
                            )
    pltt.camera.azimuth += dtheta#360/NR_FRAMES

    pltt.camera.focal_point = (0,0,0)


    pltt.show(auto_close = False)
    pltt.add_text(r"$t=%.3f$ Myr"%time,position='upper_left',
                        font_size=18,
                        color="white")
    if args.screenshot: pltt.screenshot(filename = path+"/plot_%05d.png"%i) #, window_size=[ws,ws])
    
    pltt.clear_actors()
    return
#pltt.reset_camera()
#        hello= plotter.screenshot()#("hello_%05d"%i)
#pltt.show()
    
snaplist=range(args.n1, args.n2+1,args.step)                   

if args.cpus is not None:
	with Pool(processes=args.cpus)as p:
		results=p.starmap_async(make_a_plot, zip(snaplist, repeat(args)) )
		results.get()###NEEDED TO OBTAIN THE RESULTS IN ASYNC
		p.close()
else:
	for i in snaplist:
		make_a_plot(i,args)
