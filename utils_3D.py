import pickle as pk
import os
import numpy as np
import pyvista as pv
from utils.utils import flatten
import operator as opp
verbose = True
if not verbose:
    def print(*Args, **kwargs):
        return
def clip_in(a,b):
    """Find which elements od a are in the range of values defined by b[0] and b[1]

    Args:
        a (np.ndarray): array to clip
        b (array(2)): array of values to use to perform the clip
    Returns:
        array: boolean array containg informations about which values to be kept and which not.
    """
    return np.logical_and(a > b[0], a < b[1] )
def clip_out(a,b):
    return np.logical_not(clip_in(a,b))

op ={
     ">"  : opp.gt,
     ">=" : opp.ge,
     "<"  : opp.lt,
     "<=" : opp.le,
     "eq" : opp.eq,
     "in" : clip_in,
     "out": clip_out
    
}
#pv.global_theme.allow_empty_mesh = True

def load_and_save_simple(ds,
                       quant           ,#= "magnetic_field_magnitude",
                       rs_path         = "./",
                       threshold_field = ["velocity_magnitude"],#["velocity_magnitude"],
                       axes            = ["x","y","z"],
                       load_data       = False,
                       save_data       = False):
    
    rs_path   += "/"
    dirname   = "all_data/"
    #base_name = ""
    base_path = rs_path + dirname #+ base_name
    fields = [ quant ] + axes
    
    for tr in threshold_field:
            fields += [tr]
    #print(fields)
    d=[]
    ds.init_sp()
    print_info = True
    for field in fields:
        if isinstance(field, str): 
            fname = field 
        else: 
            if field is None: fname = "empty_"+fields.index(field)
            fname = "array_"+axes.index(field)
        path = base_path + "%s.pk"%fname
        if ((load_data is True and (not os.path.exists(path)) ) or (save_data is True) or (load_data is False)):
        
            
            print("Computing data for field:", fname)
            
            
            if isinstance(field, str): 
                to_dump = ds.sp[("gas" if "hydro" not in fname else "ramses", fname)].in_cgs().d
                
            else: 
                to_dump = field
                
            d += [to_dump]
            if save_data or (load_data is True and (not os.path.exists(path))):
                
                try:
                    os.mkdir(rs_path + dirname)
                except Exception as e:
                    if print_info:
                        print(e, "impossible to create the directory.")
                        print_info = False
                with open(path, "wb") as output_file:
                    pk.dump([to_dump], output_file)
                    
        elif ((load_data is True and os.path.exists(path))):
            print("Existing data found for field: %s"%fname)
            with open(path, "rb") as output_file:
                d+=pk.load( output_file)
    #print(d[4:len(d)])
    magnitude, x, y, z = d[0:4]
    thresh_data = d[4:len(d)]
        
    return [[thresh_data,magnitude], [x, y, z]]


def load_and_save_data(ds, level,
                       quant           ,#= "magnetic_field_magnitude",
                       rs_path         = "./",
                       key             = "magnetic_field",
                       threshold_field = "velocity_magnitude",
                       load_data       = True,
                       save_data       = True,
                       vectorial_data  = True):
    rs_path   += "/"
    dirname   = "regular_grid_%02d/"%level
    base_name = ""
    fields=[
            threshold_field,
            quant
            ]
    if vectorial_data:
        fields+=[ 
            "%s_x"%key,
            "%s_y"%key,
            "%s_z"%key,
            ]
    d=[]
    print_info = True
    for field in fields:
        
        path = rs_path + dirname + base_name + "%s.pk"%field
        
        if ((load_data is True and (not os.path.exists(path)) ) or (save_data is True) or (load_data is False)):
            print("Computing data for field:", field)
            
            cube = ds.get_cube(level, field, typefield = "gas" if "hydro" not in field else "ramses" , ghost = 10)
            
            to_dump = cube[0].T
            #if "velocity" in field: to_dump.to("km/s")
            if save_data or (load_data is True and (not os.path.exists(path))):
                try:
                    os.mkdir(rs_path + dirname)
                except Exception as e:
                    if print_info: print(e, "impossible to create the directory.")
                print_info = False
                with open(path, "wb") as output_file:
                    pk.dump([to_dump.in_cgs()], output_file)
            d+=[to_dump]
            del cube
            
        elif ((load_data is True and os.path.exists(path))):
            print("Existing data found for field: %s"%field)
            with open(path, "rb") as output_file:
                d+=pk.load( output_file)
    if vectorial_data:
        thresh_data , magnitude, vec_x, vec_y, vec_z = d
    else: 
        thresh_data , magnitude = d
        
    return [[thresh_data , magnitude], [vec_x, vec_y, vec_z]] if vectorial_data else [thresh_data , magnitude]

def load_and_save_grid(level,
                       grid , 
                       thresh_mag_data,
                       vec_data,
                       key,
                       quant_name      = None,
                       rs_path         = "./",
                       load_grid       = False,
                       LOG             = True,
                       save_grid       = True):
    
    rs_path   += "/"
    dirname    = "pyv_grids_%02d/"%level
    base_name  = "streamlines_" + key +"%s"%quant_name
    
    xcut, ycut, zcut = grid
    thresh_data , magnitude = thresh_mag_data
    if quant_name is None: quant_name = "[%s]"% str(magnitude.units).replace("/","\\") ##"last_value"#"[%s]"% magnitude.name
    if LOG: quant_name = quant_name 
    path       = rs_path + dirname + base_name + "_%s.pk"%quant_name
    vec_x_cut, vec_y_cut, vec_z_cut = vec_data
    print_info = True
    if ((load_grid is True and (not os.path.exists(path)) ) or (save_grid is True) or (load_grid is False)):
        
        gridcut = pv.StructuredGrid(xcut, ycut, zcut)
        print("Grid done.")
        vectors = np.empty((gridcut.n_points, 3))
        # normalise vectors data to have values close to 1
        vflats = [vec_x_cut.flatten(),vec_y_cut.flatten(),vec_z_cut.flatten()]
        #vnorm = np.sqrt(vflats[0]**2 + vflats[1]**2 + vflats[2]**2) 
        vectors[:, 0] = vec_x_cut.flatten()#/vnorm
        vectors[:, 1] = vec_y_cut.flatten()#/vnorm
        vectors[:, 2] = vec_z_cut.flatten()#/vnorm
        print("Vectors done.")
        gridcut[key]  = vectors
        #print(vectors, key)
        #print(vec_x_cut.flatten(),magnitude.flatten())
        gridcut.point_data[quant_name] =  magnitude.flatten()# if LOG else magnitude_cut).flatten()
        #gridcut.point_data["B_magcutLOG"] = 
        gridcut.point_data["thresh_data"] = thresh_data.flatten()
        print("Data assinged to the grid.")
        #print(magnitude)
        if save_grid: 
            print("Saving the grid")
            try:
                os.mkdir(rs_path+dirname)
            except Exception as e:
                if print_info:
                    print(e, "You cannot create this directory directory.")
                    print_info = False
            with open(path, "wb") as output_file:
                        pk.dump(gridcut, output_file)
            print("saved.")
    elif load_grid is True and  os.path.exists(path):
        print("pickling") 
        with open(path, "rb") as output_file:
                    gridcut = pk.load( output_file)

    return [gridcut, quant_name]

def extract_cut(width,thresh_mag_data, vec_data, bounds = None):
    
    thresh_data , magnitude = thresh_mag_data
    vec_x, vec_y, vec_z = vec_data
    dx    = width/len(thresh_data)
    n     = np.arange(1, len(thresh_data)+1)*dx

    x, y, z = np.meshgrid(n,n,n, indexing="ij")#x, y, z
    n     = len(n)

    dx=x[1,0,0]-x[0,0,0]
    test_arr=np.array(range(1, n+1))*dx-0.5*dx
    test_ind=np.array(range(1, n+1))
    #lets define a cube based on indices
    #print(bounds)

    if bounds is not None:
        boundsx = bounds[0]
        boundsy = bounds[1]
        boundsz = bounds[2]
    else:
        boundsx = boundsy = boundsz = [0,width]
    x0,wx = boundsx[0],boundsx[1]
    y0,wy = boundsy[0],boundsy[1]
    z0,wz = boundsz[0],boundsz[1]    
    #bound = [367.18743896484375, 549.7140502929688, 241.33399963378906, 510.00592041015625, 407.23834228515625, 562.1353149414062]
    
    #print(x0, y0, z0, wx, wy,wz, test_arr)
    imin  = np.searchsorted(test_arr, x0)
    imax  = np.searchsorted(test_arr, wx)
    jmin  = np.searchsorted(test_arr, y0)
    jmax  = np.searchsorted(test_arr, wy)
    kmin  = np.searchsorted(test_arr, z0)
    kmax  = np.searchsorted(test_arr, wz)

    ##cut  
    xcut, ycut, zcut = x[imin:imax, jmin:jmax, kmin:kmax], y[imin:imax, jmin:jmax, kmin:kmax],z[imin:imax, jmin:jmax, kmin:kmax]
    magnitude_cut    = magnitude[imin:imax, jmin:jmax, kmin:kmax]
    thresh_data_cut  = thresh_data[imin:imax, jmin:jmax, kmin:kmax]
    vec_x_cut        = vec_x[imin:imax, jmin:jmax, kmin:kmax]
    vec_y_cut        = vec_y[imin:imax, jmin:jmax, kmin:kmax]
    vec_z_cut        = vec_z[imin:imax, jmin:jmax, kmin:kmax]
    
    
    return [[xcut, ycut, zcut], [thresh_data_cut , magnitude_cut], [vec_x_cut, vec_y_cut, vec_z_cut]]
    ####JUST  TO TRY, DELETE DATA TO MAKE SPACE
        
def streamplot(gridcut,
               thresh_mag_data,
               key, 
               quant_name,
               from_source = True,   
               LOG         = True,
               threshold   = 5e5,
               num_points  = 100,
               cmap        = "rainbow", 
               clip        = False, 
               radius      = None,
               max_integration_time = None,
               tf          = None,
               cut_below   = True, 
               clim        = None, 
               selected_points = None,
               pbounds      = None,
               plotter      = None,
               ):
    stream_opacity = tf is not None 
    if LOG: gridcut.point_data[quant_name] = np.log10(gridcut.point_data[quant_name])
    thresh_data_cut , magnitude_cut  = thresh_mag_data
    
    if plotter is None: plotter = pv.Plotter(off_screen = False)
    
    cond = thresh_data_cut > threshold if cut_below else thresh_data_cut < threshold 
    if clim is None:
        clim=(np.log10(magnitude_cut[cond]).min(),
            np.log10(magnitude_cut[cond]).max()) if LOG else (magnitude_cut[cond].min(),
                                                                magnitude_cut[cond].max())
    
    
    if from_source: 
        
        ## Extract the corresponding coordinates for the selected points
        if selected_points is None:
            flattened_magnitude_cut = magnitude_cut.flatten()

            ic_min,_= findclim(tf, 1, reverse = False)
            ic_max,_= findclim(tf, 1, reverse = True)
            dC      = (clim[1] - clim[0])/len(tf)
            B       = np.arange(1, len(tf)+1)
            Brange  = clim[0]+dC*(B)

            
            min_magnitude_field = 10**Brange[ic_min] if LOG else Brange[ic_min] #1e-10#10**clim[0]
            max_magnitude_field = 10**Brange[ic_max] if LOG else Brange[ic_max] #5e-10#10**clim[1]
            thresh_data_cut_f = thresh_data_cut.flatten()
            cond = thresh_data_cut_f > threshold if cut_below else  thresh_data_cut_f < threshold
            valid_indices = np.where((flattened_magnitude_cut > min_magnitude_field) &
                                    (flattened_magnitude_cut < max_magnitude_field) &
                                    (cond))[0]
            # Randomly select num_points from the valid indices
            selected_indices = np.random.choice(valid_indices, size = num_points, replace = True) ### WAS FALSE WHEN IT WORKED
            selected_points = (gridcut).points[selected_indices]
            #print(selected_points)
            # Create a PyVista PolyData object for the selected seed points
            
        seed_points_polydata = pv.PolyData(selected_points)

        # Create a PyVista PolyData object for the random seed points
        #seed_points_polydata = pv.PolyData(random_points)
        stream = gridcut.streamlines_from_source(seed_points_polydata,                                  
                                                vectors             = key,
                                                surface_streamlines = False,
                                                progress_bar        = True,
                                                max_error           = 1.0e-16,
                                                interpolator_type   = "point",
                                                max_steps           = 1000000,
                                                max_time            = max_integration_time,  
                                            )
    else:
        
        stream = gridcut.streamlines( key, n_points = num_points )

    if radius is not None: stream = stream.tube( radius = radius )

    if clip:
        stream = stream.clip_scalar(scalars = "thresh_data",
                                    invert  = not cut_below,
                                    value   = threshold)
        
    label =  r"Log10 " + quant_name if LOG else quant_name

    
    plotter.add_mesh( stream ,
                      scalars         = quant_name,
                      clim            = (clim[0],clim[1]),
                      opacity         = tf if stream_opacity else 1.,                            
                      cmap            = cmap, 
                      scalar_bar_args = {"color":"white", "title": label}
                           )
    
    plotter.show_bounds(
                        bounds    = pbounds if pbounds is None else flatten(pbounds),
                        color     = "white",
                        bold      = False,
                        location  = 'origin',                       
                        ticks     = 'both',                       
                        n_xlabels = 4,                        
                        n_ylabels = 4,                        
                        n_zlabels = 4,        
                        )
    plotter.show_axes()
    plotter.set_background("black")
    return plotter
def setup_3Ddata(
                mags,
                coords,
                scale_axes= [None, None,None],
                scaleq    = 1.,
                threshold = 0.,
                operation = [">"],
                clim      = None,
                centered= True,
                box=0,
                split=None,
                rescale=-1
                ):
    threshold_field, quant =  mags
    operation=operation*len(threshold_field)
    #print(threshold_field, operation, threshold)

    x, y, z = coords
    
    if scale_axes is not None:
        for i, scale_i in enumerate(scale_axes):
            if scale_i is not None: coords[i] /= scale_i
    if centered: 
        x-=0.5*box
        y-=0.5*box
        z-=0.5*box
    
    if scaleq is not None: quant /= scaleq
    
    if (threshold is not None) and (threshold_field is not None): 
        #print(threshold, quant)
        first = True
        count = 0
        for thres, oper in zip(threshold, operation):
            #if first: threshold_field_c = threshold_field
            trf = threshold_field[count]
            trf, quantities = clip_arrays([quant, x, y, z]+threshold_field, trf, thres, operation = oper)
            quant, x, y, z  = quantities[0:4]
            threshold_field = quantities[4:len(quantities)]
        
            count+=1
            
    if clim is not None:
        #print("ocjcijoi",x,y,z, threshold_field, quant)
        cond_min, cond_max = quant > clim[0],quant < clim[1]

        x              =              x[np.logical_and(cond_min, cond_max)]
        y              =              y[np.logical_and(cond_min, cond_max)]
        z              =              z[np.logical_and(cond_min, cond_max)]
        if threshold_field is not None: 
            if isinstance(threshold_field, list):
                for ii,_ in enumerate(threshold_field):
                    threshold_field[ii]=_[np.logical_and(cond_min, cond_max)]
        quant          =          quant[np.logical_and(cond_min, cond_max)]
    if split is not None:
        
        if split[0]=="x" : cutter = x.copy()
        if split[0]=="y" : cutter = y.copy()
        if split[0]=="z" : cutter = z.copy()
   
        cond = op[split[2]](cutter,split[1]) #_min, cond_max 
        x              =              x[cond]
        y              =              y[cond]
        z              =              z[cond]
        if threshold_field is not None: 
            if isinstance(threshold_field, list):
                for ii,_ in enumerate(threshold_field):
                    threshold_field[ii]=_[np.logical_and(cond,1)]
            else:
                
                threshold_field=threshold_field[np.logical_and(cond,1)]
        quant          =          quant[cond]
        del cutter
    if rescale>-1:
        for jjj in range(rescale):
            x              =              x[1::2]
            y              =              y[1::2]
            z              =              z[1::2]
            if threshold_field is not None: 
                if isinstance(threshold_field, list):
                    for ii,_ in enumerate(threshold_field):
                        threshold_field[ii]=_threshold_field[1::2]
            quant          =          quant[1::2]
    return [[threshold_field, quant], [x,y,z]]

def clip_arrays(quantities, threshold_array, threshold, operation):
    
    if not isinstance(quantities, list): quantities=list(quantities)
        
    cond =  (op[operation])(threshold_array, threshold)
                
    for i,q in enumerate(quantities):
        quantities[i] = q[cond]
    threshold_array = threshold_array[cond]
    
    return threshold_array, quantities
    
def bins3d(x,y,z, bins = 128, density=False):
    t =  np.zeros_like((x,y,z)).T
    t[:, 0]= x
    t[:, 1]= y
    t[:, 2]= z    
    H, edges = np.histogramdd(t, bins=bins, density=density)
    for idx, edge in enumerate(edges): 
        edges[idx] = 0.5*(edge[1:]+edge[:-1])
    return H, edges

def plot3D(coords, 
           qname, 
           quant,
           pbounds    = None, 
           point_size = 0.8,
           render_points_as_spheres = False, 
           gaussian = False, 
           log = True, 
           diffuse  = None, 
           tf       = None, 
           cmap     = "rainbow", 
           plotter  = None, 
           labels   = ["X axis","Y axis","Z axis"]
           ):
    if plotter is None: plotter = pv.Plotter(off_screen = True)
    style = "points" if not gaussian else "points_gaussian"
    if tf is None: tf = 1.
    if diffuse is None: diffuse = 1 if render_points_as_spheres else 0.03/5. #0.07#0.07
    
    x,y,z = coords 
    
    points = np.column_stack((x, y, z))
    point_cloud = pv.PolyData(points)
    #print(x,y,z, quant)
    point_cloud[qname] = quant
    # Create a PyVista plotter


    print("Building plot ")
    render_points_as_spheres=False
    print(render_points_as_spheres, gaussian, style)
    plotter.add_points(
                    point_cloud,
                    scalars=qname,
                    point_size = point_size,
                    cmap      = cmap,
                    style     = style,
                    render_points_as_spheres = render_points_as_spheres,
                    emissive  = gaussian ,
                    log_scale = log,
                    opacity  = tf,
                    diffuse   = diffuse,

                    scalar_bar_args={"color":"white"}
                    )

    plotter.set_background("black")
    plotter.show_bounds(    bounds=pbounds if pbounds is None else flatten(pbounds),
                            color="white",
                            bold=False,
                            location='outer',                       
                            ticks='both',                       
                            n_xlabels=4,                        
                            n_ylabels=4,                        
                            n_zlabels=4,                        
                            xtitle=labels[0],                       
                            ytitle=labels[1],                      
                            ztitle=labels[2],    
                            font_size = 20,
                            
                            ##font_size=20                 
                            )


    return plotter













def findclim(gmed,
             treshold,
             reverse = False):
    """This function find the index of the first element of an array exceeding a certain threshold. 

    Args:
        gmed     (array):
                        Array to examine.
        treshold (float, int): 
                        Threshold valeìue. 
        reverse (bool, optional):Start from the end of the array going backwards. Defaults to False.

    Raises:
        Exception: General errors.

    Returns:
        (int, float): Index and values of the array at the exceeding point.
    """
                
    cell_center=np.arange(1,len(gmed)+1 )
    if reverse:
        gmed=np.flip(gmed)
        cell_center=np.flip(cell_center)
    index=0
    
    for i in gmed:
        ciplux = True
        if i >= treshold and np.isnan(i) == False and np.isinf(i) == False: 
            #print(i, "eccolo!!", mean_value, reverse_shock_pos_index)
            ciplux=False
            break
        if ciplux: index += 1
    try:
        out = len(gmed) - index - 1  if reverse else index
        if out < 1: out = 1 
        if out > len(gmed)-1: out = len(gmed)-1  
        return out,  gmed[index] #cell_center[index]
    except:
        raise Exception("Cant build the colorbar")### se non lo trova arriva all'inizio, quindi metto semplicemnte 0.
def set_diffusive(quant):
    
    l = len(quant)
    
    fac = 10
    n = np.log10(l)**2
    print(n)
    return 3/(fac*n) #if n>5 else 0.3

def set_size(quant, size, min_size):
    n=np.log10(len(quant))
    slope = -1
    y     = size+slope*np.log10(n)
    return np.array([y,min_size]).max()
