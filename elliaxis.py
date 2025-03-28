#Moment of inertia tensor (center_of_mass must not have any units attached to it!)
import yt
import numpy as np 
import analysis_playground.analysis_tools as an
import analysis_playground.constants as cc

import argparse
from multiprocessing import Pool


def clean_arrays(arrays, condition):
    for n,array in enumerate(arrays):
        arrays[n] = array[condition]
    return arrays

# Calculate the moment of inertia tensor
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
parser.add_argument("-n2", 	help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None)
parser.add_argument("-out","--output",	  help ="Set the output file name.", default="axes_cm_mtot.txt",type = str)
parser.add_argument("-vlim",help ="Minimum velocity.", default = None, type=float,      nargs=1)
parser.add_argument("-inklim",help ="Minimum velocity.", default = None, type=float,      nargs=1)
parser.add_argument("-templim",help ="Minimum velocity.", default = None, type=float,      nargs=1)
parser.add_argument("-maglim",help ="Minimum velocity.", default = None, type=float,      nargs=1)
parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)
parser.add_argument("-cpus", 	help="step in the data. Use only if n2 is used.", type=int, default=None)
#
args = parser.parse_args()
if args.n2 is None: args.n2=int(args.n1)
if  args.n1 > args.n2: 
    n2tmp	= int(args.n2)
    args.n2 = int(args.n1)
    args.n1 = n2tmp
snaplist=range(args.n1, args.n2+1,args.step)                   

if args.cpus is not None:
    print("importin multirpocess")
    from multiprocessing import Pool
    from itertools import repeat

def cmp_center_of_mass(xx, yy, zz, mm):
    MM = np.sum(mm)
    xm = np.sum(mm*xx)/MM
    ym = np.sum(mm*yy)/MM
    zm = np.sum(mm*zz)/MM
    return [xm, ym, zm], MM

def compute_inertia_tensor(x, y, z ,mass, center_of_mass):
    # Extract data fields
    #cloud_patch_mass = region["gas", "cloud_patch_mass"].v
    dx = x - center_of_mass[0]
    dy = y - center_of_mass[1]
    dz = z - center_of_mass[2]
    
    # Calculate inertia tensor components
    I_xx = np.sum(mass * (dy*2 + dz*2))
    I_yy = np.sum(mass * (dx*2 + dz*2))
    I_zz = np.sum(mass * (dx*2 + dy*2))
    
    I_xy = -np.sum(mass * dx * dy)
    I_xz = -np.sum(mass * dx * dz)
    I_yz = -np.sum(mass * dy * dz)

    inertia_tensor = np.array([[I_xx, I_xy, I_xz],
                               [I_xy, I_yy, I_yz],
                               [I_xz, I_yz, I_zz]])

    return inertia_tensor
    
def make_axes(i, args):
    print("Initialise data.")
    ds = an.analyse(outnumb = i)
    ds.init_sp()
    time = ds.time.to("kyr").v
    sp = ds.sp
    x = sp["gas","x"].v/cc.pc
    y = sp["gas","y"].v/cc.pc
    z = sp["gas","z"].v/cc.pc
    mass = sp["gas","mass"].v/cc.M_sun
    arrays = [x, y ,z , mass]
    print("Cleaning data.")
    if args.vlim is not None:   
        arrays = clean_arrays(arrays, sp["gas","velocity_magnitude"].in_cgs() > args.vlim )
    if args.templim is not None: 
        cond= sp["gas","temperature"] > args.templim
        x    = x   [cond]
        y    = y   [cond]
        z    = z   [cond]
        mass = mass[cond]
        #arrays = clean_arrays(arrays, sp["gas","temperature"]              > args.templim )
    print(max(x),max(y),max(z),min(x),min(y),min(z))
    print("Arrays are clean, computing mass.")  
    center_of_mass, mtot = cmp_center_of_mass(x, y, z, mass)
    print("center of mass",center_of_mass, mtot)
    print("computing initertia tensor.")
    I_tensor = compute_inertia_tensor(x, y, z ,mass, center_of_mass)
    # Diagonalize the moment of inertia tensor
    eigenvalues, eigenvectors = np.linalg.eigh(I_tensor)
    # Calculate the semi-axes lengths
    a_squared = (5.0 / (2.0 * mtot)) * (eigenvalues[1] + eigenvalues[2] - eigenvalues[0])
    b_squared = (5.0 / (2.0 * mtot)) * (eigenvalues[0] + eigenvalues[2] - eigenvalues[1])
    c_squared = (5.0 / (2.0 * mtot)) * (eigenvalues[0] + eigenvalues[1] - eigenvalues[2])
    print(a_squared,b_squared,c_squared)
    axes_lengths = np.sqrt([a_squared, b_squared, c_squared])
    with open(args.output, mode="a") as f:
        txt = " % i %e %e %e %e %e \n"%(i, time, axes_lengths[0],axes_lengths[1],axes_lengths[2], mtot)
        f.write(txt)
    print("done")
    return 

if args.cpus is not None:
	with Pool(processes=args.cpus)as p:
		results=p.starmap_async(make_axes, zip(snaplist, repeat(args)) )
		results.get()###NEEDED TO OBTAIN THE RESULTS IN ASYNC
		p.close()
else:
	for i in snaplist:
		make_axes(i,args)

#if args.inklim is not None: 
#    for i in range(0,2+1): 
#        add_pscalar(i)
#    arrays = clean_arrays(arrays, sp["gas","hydro_scalar_00"]          > args.inklim )
#    
#
#    
#if args.maglim is not None: 
#    arrays = clean_arrays(arrays, sp["gas","magnetic_field_magnitude"].in_cgs() > args.maglim )