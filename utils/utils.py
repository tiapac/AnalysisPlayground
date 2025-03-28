import numpy  as np
# from CuloPy.reusable.get_info import get_info  
from scipy.io import FortranFile
FortranFile.read_vector = FortranFile.read_record
import constants as cc
import time as Time
import builtins

"""
    Here there are some usefull funtions.
"""
import inspect 



grey = "\x1b[38;20m"
yellow = "\x1b[33;20m"
red = "\x1b[31;20m"
bold_red = "\x1b[31;1m"
purple = "\x1b[35;20m"
green = "\x1b[32;20m"
blue = "\x1b[34;20m"
cyan = "\x1b[36;20m"
reset = "\x1b[0m"

def printg(args):
        print(green + args + reset)

def printgrey(args):
        print(grey + args + reset)

def printyellow(args):
        print(yellow + args + reset)

def printr(args):
        print(red + args + reset)

def printboldr(args):
        print(bold_red + args + reset)

def printpurple(args):
        print(purple + args + reset)

def printblue(args):
        print(blue + args + reset)

def printcyan(args):
        print(cyan + args + reset)


       #def __call__(self, *args, **kwargs):
       


def secure_eval(*args, **kwargs):
       for arg in args:
              for noword in ("rm", "os", "import"): 
                if noword in arg: raise Exception("'%s' is not allowed in eval."%noword)
       return eval (*args, **kwargs)





def checkstack():
        res = inspect.stack()[1]#[3]
        ini = ""#("----" if res !="makeplot" else "")
        return ini+f"called by «{res}»\t"

class stopwatch:
    def __init__(self, name = None):
        self.name = None
        self.initime =  Time.time()
        self.lasttime = float(self.initime)
        pass
    def __call__(self, msg = ""):
        right_now = Time.time()
        print("--- %s %.2e seconds ---" %(msg,(right_now - self.lasttime)))
        self.lasttime = right_now
        return
    
    def __repr__(self):
        pass


def oneD_avg_weighted_bins(x, y, w = None, bins = 128, maxrange = None):
        if maxrange is None: maxrange = [0,np.max(x)]
        
        with np.errstate(divide='ignore',invalid='ignore'):
                    gmed, edges   = np.histogram(x, bins = bins, range = maxrange, weights = y if w is None else y * w)
                    ncelle, edges = np.histogram(x, bins = bins, range = maxrange, weights = w)
                    cell_center = 0.5*(edges[1:]+edges[:-1])
                    gmed /= ncelle
        return  cell_center,gmed

def find_axes(axis, select_from_array = None):
        "This given  axis == 'x','y' or 'z' in 'xyz', this function returns the other two in alphabetical order."
        all_axes = ["x","y","z"]
        
        all_axes.remove(axis)
        if select_from_array is None: 
                
                return [min(all_axes), max(all_axes)]
        else:   
                dictionary = {
                        "x":select_from_array[0],
                        "y":select_from_array[1],
                        "z":select_from_array[2]
                        }
                return dictionary[min(all_axes)],dictionary[max(all_axes)]
def weighted_avg_and_std( values, weights):
        """
        Return the weighted average and standard deviation.

        values, weights -- NumPy ndarrays with the same shape.
        """
        average = np.average(values, weights=weights)
        # Fast and numerically precise:
        variance = np.average((values-average)**2, weights=weights)
        return (average, np.sqrt(variance))

def add_xray_stuff_AREPO(ds):
        ftype="gas"
        ptype="PartType0"
        unit_system = ds.unit_system
        #for i in ds.derived_field_list: print(i)
        from yt.utilities.chemical_formulas import ChemicalFormula
        from yt.utilities.physical_ratios import _primordial_mass_fraction

        A_H = ChemicalFormula("H").weight
        m_u = ds.units.physical_constants.amu_cgs
        X_H = _primordial_mass_fraction["H"]

        def _h_number_density(field, data):
                                return data[ptype, "density"] * X_H / (A_H * m_u)

        def _h_number_density(field, data):
                                return (
                                    data[ptype, "density"]*
                                    X_H #* data["gas", "H_fraction"]
                                    / (A_H * m_u)
                                )
        ds.add_field(
                (ptype, "H_number_density"),
                sampling_type="particle",
                function=_h_number_density,
                units=ds.unit_system["number_density"],
                        )

        def _el_number_density(field, data):
                return (
                        data[ptype, "ElectronAbundance"] * data[ptype, "H_number_density"]
                )

        ds.add_field(
                (ptype, "El_number_density"),
                sampling_type="particle",
                function=_el_number_density,
                units=unit_system["number_density"],
                )

        def _emission_measure(field, data):
                #print("hello hellooo")
                dV = data[ptype, "mass"] / data[ptype, "density"]
                nenhdV = data[ptype, "H_nuclei_density"] * dV
                nenhdV *= data[ptype, "El_number_density"]
                return nenhdV
        ds.add_field(
                (ftype, "emission_measure"),
                sampling_type="local",
                function=_emission_measure,
                units=unit_system["number_density"],
            )
        
        return
    

def add_some_fields(ds):
        fields = [("gas", "density"),
                ("gas", "temperature"),
                ("gas", "velocity_x"),
                ("gas", "velocity_y"),
                ("gas", "velocity_z"),
                #("gravity", "Potential")
                ]
        grad_fields = ds.add_gradient_fields(fields)

        # add velocity divergence. 
        def __myvdiv(field, data):
                return (
                (data['gas', 'velocity_x_gradient_x']+data['gas', 'velocity_y_gradient_y']+data['gas', 'velocity_z_gradient_z'])
                )
        ds.add_field(
                name=("gas", "mydiv"),
                function=__myvdiv,
                sampling_type="local",
                display_name=r"$\nabla\cdot v$"
                
                )
        return



def flatten(xss):
    return [x for xs in xss for x in xs]


def set_axes_equal( ax):
        """
        Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc.

        Input
        ax: a matplotlib axis, e.g., as output from plt.gca().
        """

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5*max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
        return ax



def weighted_profile(q, x, w, bins = 512, range = None):
        with np.errstate(divide='ignore',invalid='ignore'):
            qmed,   edges = np.histogram(x, bins=bins, range=[0,np.max(x)], weights=q*w)
            ncelle, edges = np.histogram(x, bins=bins, range=[0,np.max(x)], weights=w)
            cell_center = 0.5*(edges[1:]+edges[:-1])
            qmed /= ncelle
        return qmed, cell_center
    
    
class code_scales:
    unit     = 1.0
    scale_l  = 1.0
    scale_d  = 1.0
    scale_t  = 1.0
    def __init__(self, info = None):
        if not None:
            self.scale_l = info["unit_l"]
            self.scale_d = info["unit_d"]
            self.scale_t = info["unit_t"]
        self.scale_m      = self.scale_d * self.scale_l**3   
        self.scale_v      = self.scale_l / self.scale_t 
        self.scale_nH =  0.70651e0  / cc.mH * self.scale_d
        self.scale_b = (self.scale_t /(np.sqrt(4.* np.pi* self.scale_d) *self.scale_l))**(-1) 
        pass
    def __getitem__(self, key):
        return getattr(self, key)

# def load_map(dir, proj, kind, num ):
    
#     """Load a fortran binary 2D map/movie RAMSES output"""
#     # define map path
#     map_file = "%s/movie%d/%s_%05d.map" % (dir, proj, kind, num)
    
#     info    = get_info("%s/movie%d/info_%05d.txt" % (dir,proj, num))
        
#     # read image data
#     with  FortranFile(map_file) as ffile:
#         dat0 = ffile.read_vector('d')
#         [frame_nx, frame_ny] = ffile.read_vector("i")
#         dat = ffile.read_vector('f4').reshape(frame_nx,frame_ny )#, dtype=np.float64
    
    
#     return [dat, dat0, [frame_nx, frame_ny], info]

def find_shock_reverse(gmed, cell_center, treshold, mean_value, invert = False, smaller = False):###treshold in percentual
        
    reverese_gmed       = gmed if invert else np.flip(gmed)
    reverese_cell_center= cell_center if invert else np.flip(cell_center)
    reverse_shock_pos_index = 0
    treshold_toll_reverse   = 1.+treshold/100.
    for i in reverese_gmed:
        ciplux=True
        tresholdc = treshold_toll_reverse*mean_value
        cond = (i>=tresholdc) if not smaller else (i<=tresholdc)
        #print(i, tresholdc)
        if cond and np.isnan(i)==False and np.isinf(i)==False: 
            ciplux=False
            break
        if ciplux: reverse_shock_pos_index+=1
    try:
        return reverese_cell_center[reverse_shock_pos_index]
    except:
        return 0.

        # reverese_gmed = gmed if invert else np.flip(gmed)
        # reverese_cell_center = cell_center if invert else np.flip(cell_center)
        # treshold_toll_reverse = 1. + treshold / 100.

        # valid_indices = np.where((reverese_gmed >= treshold_toll_reverse * mean_value) & 
        #                                                 ~np.isnan(reverese_gmed) & 
        #                                                 ~np.isinf(reverese_gmed))[0]

        # if valid_indices.size > 0:
        #         return reverese_cell_center[valid_indices[0]]
        # else:
        #         return 0.
