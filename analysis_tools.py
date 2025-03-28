import yt 
import constants as costs
import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.axes_grid1 import AxesGrid
from other import *
import os 
from matplotlib.patches import Rectangle
import paper_settings as pg
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter, NullFormatter
import matplotlib.path as mplPath
from scipy.ndimage.filters import gaussian_filter
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable


paper = pg.PaperSettings()

myunits = paper.units
plt.rcParams.update({
			"text.usetex": True, **paper.mplfigfontstyle
			#"font.family": "Helvetica"
			})
try:
    from pyvista import Plotter
    from pyvista import PolyData
    import utils_3D
    import transfer_functions as tff


except: 
    nopyvista = True
    pass
try:
    from mylibraries.addedfields import initfield
    initfield_exitsts = True
    initfield(MHD=True)
except Exception as e :
    print(e, "File 'mylibraries.addedfields' not found.") 
    initfield_exitsts = False

import unyt
import utils.utils as utils
#units_override = {"velocity_unit": (1, "km/s"), "length_unit":(1, "parsec"),"mass_unit":(1, "Msun") }
scale_l = 3.0856776e18 
factG_in_cgs = 6.67259e-8
scale_m= 1.9891000e33
scale_d = scale_m/scale_l**3
#scale_t converts time from user units into seconds
scale_t = 1./np.sqrt(factG_in_cgs*scale_d)
scale_v = scale_l / scale_t
units_override = {"time_unit": (scale_t, "s"), "length_unit":(scale_l, "cm"),"mass_unit":(scale_m, "g") }

yt.set_log_level(50)
  

class multi:
    def __init__(self, path="./"):
        self.ts = yt.load(path+"output_?????")#/info_?????.txt")
        self.nsnaps=len(self.ts)

    def data(self):
        return self.ts, self
    
    def particle_trajectories(self, plot=False, marker=r"$\odot$", savefig=True, fig=None,ax=None, weight=("mass", "Msun",100), cunits="parsec"):
        x,y,z=[],[],[]
        first=True
        #print(enumerate([analyse(outpath=ds) for ds in self.ts ]))
        #sys.exit()
        i=1
        #print(self.ts)
        for  ds in [analyse(outpath=ds) for ds in self.ts ]:  #self.ts:#[analyse(outpath=ds) for ds in self.ts ]:
            xx,yy,zz, qq = ds.particle_data(quant="particle_index", cunits=cunits)
            cl = ds.particle_data(quant="particle_cluster_p", only_data=True)
            #if not(isinstance(weight,(int, float, tuple))) and not (weight is None): raise Exception("weight cannot be accepted.") 
            if weight is not None:
                if isinstance(weight,int) or isinstance(weight,float):
                    numweight=True
                    pass
                else:
                    ww = ds.particle_data(quant="particle_"+weight[0],only_data=True,qunits=weight[1])
                    if isinstance(weight[2],int) or isinstance(weight[2],float): 
                        ww=ww*weight[2]
                    elif hasattr(weight[2], '__call__'):
                        ww=weight[2](ww)
            if first: 
                npart=len(xx)
                x=np.zeros((self.nsnaps, npart))
                y=np.zeros((self.nsnaps, npart))
                z=np.zeros((self.nsnaps, npart))
                q=np.zeros((self.nsnaps, npart))
                cluster=np.zeros((self.nsnaps, npart), dtype=int)
                if weight is not None: w=np.zeros((self.nsnaps, npart), dtype=int)
                first=False
                
            p=np.argsort(qq)

            x[i-1,:]=xx[p]
            y[i-1,:]=yy[p]
            z[i-1,:]=zz[p]
            cluster[i-1,:]=cl[p]
            q[i-1,:]=qq[p]
            #print(x[i-1,:])
            if (weight is not None) and numweight is False: w[i-1,:]=ww[p]
            i+=1
        if plot: 
            from mylibraries.funclist import colori
            #if fig==None: 
            fig=plt.figure(figsize=(8,8))
            #if ax==None: 
            ax=fig.add_subplot(111,projection="3d")
            
            for i in range(0,npart): 
                color=colori(1, seed=cluster[0,i]) 
                ax.plot(x[:,i],y[:,i],z[:,i], c=color[0], alpha=1)
                if weight!=None: 
                    ax.scatter(x[self.nsnaps-1,i],
                               y[self.nsnaps-1,i],
                               z[self.nsnaps-1,i],
                               marker=marker,
                               s=w[self.nsnaps-1,i],
                               c=color[0], alpha=1)
        
        return x,y,z,cluster, q
        
    def __iter__(self):
        return [analyse(outpath=ds) for ds in self.ts ]
    def __len__(self):
        return len([i for i in self.ts])
    def __repr__(self) -> str:
        return "%r"%[i for i in self.ts]
    


class analyse:
    """
        This class simplify analysis of simulation snapshots. Currently tested only on RAMSES and partially on AREPO.
        outpath      :     str|yt dataset, path to the simulations output or a yt dataset. It can be the complit path or the path to the directory
                                       where the data is stored.
        base_filename: str,            base of the output. defaults to output_, which is the case for RAMSES
        fmt          : str,            numbering format of  the ouputs, defaults to "%05d".
        extention    : str,            extention to add to the output name. Defaults to "", which indicates a directory.
        outnumb      : int,            number of the simulation output
            example: In the AREPO case a dataset can be loaded in two ways: ds = analyse(filepath)
                                        or as 
                                        analyse("path_to_dir",outnumb=<int>, base_filename = "snapZ_5.000000000000001_", fmt = "%03d", extention=".hdf5", smoothing_factor=3)

        
        time_series  : None,           this is a dummy variable for now
        mhd          : bool,           if the output contains mhd data, it should be se to True.
        select_data  : tuple ("in"|"out", (fieldtype, field), vmin, vmax): allowes you to extract only the data in the range 
                                                                         defined by these parameters.
        more fields  : bool , add some pre-defined user defined fields in utils.
        smoothing_factor: arepo cells are threated as sph particles to which is associated a smoothing lenght computed
                          as the equivalent spherical radius of the cell multiplied by "smoothing_factor".
                          Defaults to 2 in yt, here is set to 3.
                          More at https://yt-project.org/doc/examining/loading_data.html#loading-sph-data
    """

    def __init__(self, outpath = "./",
                 base_filename = "output_",
                 fmt           = "%05d", 
                 extention     = "",
                 outnumb:int = -1,
                 #time_series:bool = False, 
                 mhd:bool = False,
                 poisson=False,
                 select_data = None,
                 more_fields = False,
                 cached = True,
                 #smoothing_factor = 3. # used only for arepo when data are being loaded in the class
                 **kwargs
                 ):
        self.cached = cached
        self.outnumb=outnumb
        self.pathinfo = {"outpath":outpath, 
                         "outnumb":outnumb,
                        "base_filename":base_filename,
                        "fmt":fmt, 
                        "extention":extention}
        # mhd settings
        self.mhd     = mhd
        self.poisson = poisson
        # initialise additional fields from function init field
        if initfield_exitsts:       
            initfield(MHD = self.mhd, POISSON = self.poisson)  
        
        # load dataset 
        if  isinstance(outpath,str): 
            # first we try to load from the outpath string, in case this directly point to the dataset.
            # if it doesn't work, use also base_filename, fmt and number of the output.
            # if cache:
            #         try:
            #             self.ds = pickle.load(open(outpath, "rb"))
            #         except:
            #             self._load_ds()
            if False:
                print("Checking cached data")
                filename = (self.pathinfo["outpath"] + "/"
                    + self.pathinfo["base_filename"]
                    + self.pathinfo["fmt"]%self.pathinfo["outnumb"]
                    + self.pathinfo["extention"])
                filenamepk = filename + "/ds.pk"
                try:
                    
                    self.ds = pickle.load(open(filenamepk, "rb"))
                    print("Loaded cached data from %s"%filenamepk)
                except Exception as e:
                    print(e)
                    print(f"No cached data found. Loading data from {filename}")
                    self._load_ds(outpath)
                    pickle.dump(self.ds, open(filenamepk, "wb"))
                    print("Data saved in %s"%filenamepk)
            else:
                self._load_ds(outpath)    
        else:
            # in case a datases is already supplied 
            self.ds = outpath
            outpath = None
        
        self.directory = self.ds.directory
        self.filename  = self.ds.filename
        self.output_identifier = os.path.splitext(self.filename)
        #print(self.ds.smoothing_factor)
        # save the dataset type (rames, arepo ...) ref = int(np.prod(ds.ref_factors[0:max_level]))
        self.dataset_type = self.ds.dataset_type
        self.arepo = False
        if self.dataset_type == "arepo_hdf5": self.arepo = True
        # define gradients this have been tested only FOR RAMSES 
        if more_fields: 
            utils.add_some_fields(self.ds)
        
        self.time        = self.ds.current_time                             # current time of the simulation
        self.fields      = self.ds.derived_field_list                       # save list of available fields
        #self.ds.define_unit("Be", (1.e51, "erg"),  tex_repr=r"Be")         # add Be = 1e51 ergs as units to ds
        self.width       = self.ds.domain_width.to("parsec")                # width of the simulation box, in parsec
        
        #####################################################################################################
        # might not work with non-chartesian grid codes
        if not self.arepo:
            self.maxlevel    = self.ds.index.max_level                          # maximum level of refinment
            self.base_level  = int(np.log2(self.ds.domain_dimensions[0]))       # minimum refinement level
            self.dx_max      = (self.width.v)[0]/self.ds.domain_dimensions[0]   # minimum resolution achievable
            self.dx_min      = self.dx_max / 2**self.maxlevel                   # maximum           ""
        #####################################################################################################

        # initialise variables and set them based on inputs
        
        self.select_data = select_data                                          
        self.sp          = None
        self.data_restr  = None 

        # useful to make computation only once in some cases
        self.thermal_eng = None
        self.kinetic_eng = None
        self.tot_eng     = None
        self.radial_mom  = None
        self.x = None
        self.y = None
        self.z = None

    def _load_ds(self, name):
        try:     
                self.ds = yt.load(+ self.pathinfo["outpath"], 
                              #smoothing_factor = smoothing_factor,
                              default_species_fields='ionized'
                              )  
        except:
            
            filename = (self.pathinfo["outpath"] + "/"
                    + self.pathinfo["base_filename"]
                    + self.pathinfo["fmt"]%self.pathinfo["outnumb"]
                    + self.pathinfo["extention"])    
            self.ds = yt.load(filename, 
                            #smoothing_factor = smoothing_factor,
                            #units_override=units_override
                            default_species_fields='ionized',
                            
                            )   
    def init_xray_fields(self, 
                        kev_range    = [0.05,11.0],
                        max_density  = (5.0e-25, "g/cm**3"),
                        t_range      = None,#[ 0.025,64]
                        Zmet         = "auto",
                        spectrum_bins= 1000, # or just a list of two values  ie [0.05,1]
                        nbands       = 5,
                        redshift = 1e-5,
                        ):
        import pyxsim

        """
        This function initiate x-ray fields from pyxsim that can then be plotted or used to analyse data.

        Args:
            kev_range (list, optional)   : range of energies to consider for the spectrum. Defaults to [0.05,11.0].
            max_density (tuple, optional): Denisties over this threshold are not included. Defaults to (5.0e-25, "g/cm**3").
            t_range (_type_, optional)   : Only gas in this range of temperature is included.Units must be in Kelvin.
                                           If None, defaults to (0.025 Kev,64 Kev)kB_SI/e. 
            Zmet (Float|str, optional)   : Metallicity of the ISM, in solar units. if "auto" it tries to compute an 
                                           average from the dataset. Defaults to "auto".
            spectrum_bins (int, optional): How many bins for the emission spectrum model. Defaults to 1000.
            nbands   (int/list, optional): How many x-ray bands should be generated between between kev_range values.
                                           Defaults to 5. If <2 it will be set as 2 and the results will be just one band.
            redshift (float)              : redshift of the source. Defaults to 1e-5.                               
        Returns:
            list: list of the generated x-ray fields.
        """
        kev_range = np.array(kev_range)
        if t_range is None: t_range   = np.array([0.025,64 ])/(costs.kB_SI/costs.e_SI)
        t_range = np.array(t_range)
        # Build for arepo fields. This is very problematic and it is basically forcing some things.
        # Probably the yt-project devs should be contacted. 
        if self.dataset_type != "ramses": utils.add_xray_stuff_AREPO(self.ds)
        
        # convert temperature range to Kev
        kt_range = t_range * (costs.kB_SI / costs.e_SI)
        
        # try to compute average metallicity of the ISM
        if Zmet=="auto":
            if self.sp is None: self.init_sp()
            try:
                #self.sp.quantities.weighted_average_quantity(field, weight)
                Zmet = self.sp.quantities.weighted_average_quantity(("gas","metallicity"),
                                                                    ("gas", "mass"))
            except Exception as e: 
                # if no  metallicity field is found or given, just assume 0.
                print(e, "No info found about gas metallicity. Assuming Z = 0")
                Zmet=0.
        
        # build the thermal soruce model.
        source_model = pyxsim.CIESourceModel("cloudy",
                                             emin   = kev_range[0], emax = kev_range[1], 
                                             kT_min = kt_range[0], kT_max = kt_range[1],
                                             nbins  = spectrum_bins, 
                                             max_density=max_density,
                                             Zmet   = Zmet,
                                             
                                             binscale = "log" )
        self.thermal_sourceModel = source_model 
        # add the fields to yt. 
        if isinstance(nbands, int):
            # here define several several bands
            xray_fields =[]
            if nbands < 2: nbands = 2 
            bands = np.geomspace(kev_range[0], kev_range[1], num = nbands)
            
            for i in range(0, len(bands)-1): 
                
                xray_fields += source_model.make_source_fields   (self.ds, bands[i], bands[i+1]  )
                xray_fields += source_model.make_intensity_fields(self.ds, bands[i], bands[i+1]  ,redshift=redshift  )
        else:
            # only one band corresponding to the kev_range
            xray_fields  = source_model.make_source_fields(   self.ds, kev_range[0], kev_range[1]  )
            xray_fields += source_model.make_intensity_fields(self.ds, kev_range[0], kev_range[1], redshift=redshift  )
            #make_intensity_fields
        return xray_fields, source_model
    def emission_spectrum(self, nbins=5000, redshift = 1e-5, fig = None, ax = None, **kwargs ):
        """This function create the spectrum of a source.
        Before calling this function, a model should be created through init_xray_fields.
        Args:
            nbins    (int  , optional): Number of bins for the spectrum. Defaults to 5000.
            redshift (float, optional): Redshift of the souce. Should match the one used for the mode.
            
        
        Returns:
            list: List containing the plot axis, figure, energy bins centers and fluxes.
        """
        raise Exception("Not working for some reasons.")
        fig = plt.figure() if fig is None else fig
        ax  = fig.add_subplot(1,1,1) if ax is None else ax
        
        model = self.thermal_sourceModel
        #print(model, self.sp)
        emin  = model.emin
        emax  = model.emin 
        if self.sp is None: self.init_sp()
        pspec = model.make_spectrum(self.sp, emin, emax, nbins , redshift = redshift)        
        ax.loglog(pspec.emid, pspec.flux, **kwargs)
        return [ ax, fig, pspec.emid, pspec.flux]



###PLOTS
    def Splot(self, quant, axis = "x", typefield = "gas", typeweight = "gas", weight = None,
              log = True, units = None,proj=False,
              width = None, width_units="pc",center = "c",#[0.5,0.5,0.5],
              particles = True, p_size = ("mass", "Msun",10), marker = r"$\odot$",
              stream_v = False, stream_m = False,  cmap = None, 
              vval  = None, 
              hide_axes = False, time  = True, scale = False, grid  = False):
        """
            This function produces Slices of a simulation output.
            quant       : str  , quantity to plot
             
            weight      : str  , optional, weight used to make projections. Only used for projections
            axis        : str  , axis to slice
            typefield   : str  , type of the field of the quantity
            typeweight  : str  , type of the weight
            
            log         : bool , use log scale for colormap
            units       : str, units to use to scale the field plot and colormap 
            proj        : bool , Defaults to slice (False).
            width       : float, width of the region to plot
            width_units : str  , units to be used for the width
            center      : list[float*3], where to center the slice, axes range from 0 to 1
            
                        
            particles   : bool, plot particles positions
            p_size      : float|int|(field, units, int|float|func): give rules to scale the size of the particle markers
                        THIS ACTUALLY REQUIRES TO PATCH THE YT FUNTION
            marker      : str, select the marker for the particles positions 
            
            stream*     : bool, plot streamlines of the velocity (v) or magnetic_field(m) field 
            cmap        : str, select a different color map 
            vval        : tuple|list|np.array, select only a range of values.
            hide_axes   : bool, do not show axes of the plot
            time        : bool, show simulation time on the plot
            scale       : bool, show a scale bar in the plot
            grid        : shock the cell edges 
        """
        # set width if not supplied
        if width is None: width = float(self.width[0].d)
        
        # initialise data if a filter is applied
        if self.select_data is not None: self.init_sp()
        
        # define quantity tuple
        quant = (typefield, quant)
        
        # produce the plot with yt
        if not proj:
            s = yt.SlicePlot(self.ds, axis, quant, width=(width,width_units), center = center, data_source = self.sp if self.data_restr is not None else None  )
        else:
            if weight is not None: weight = (typeweight, weight)
            s = yt.ProjectionPlot(self.ds, axis, quant, weight_field = weight, width = (width,width_units), center = center, data_source = self.sp if self.data_restr is not None else None  )
             
        # yt does not automatically flip the plot here. We do it. 
        axses = utils.find_axes(axis)
        if axis == "y":
            s.swap_axes()
            axses.reverse()
        #    s.flip_vertical()
        
        # additional stuff
       
        if not log           : s.set_log(quant, log = False)
        if time              : s.annotate_timestamp()
        if particles         : 
            if self.ds.dataset_type=="arepo_hdf5":
                print("partilces annotation does not work really well with AREPO.\
                      \n This is because AREPO data are threated as sph data, and therefore every cell is then a particle. ")
            try:
                s.annotate_particles(width=(width,width_units), marker = marker, p_size = p_size)
            except Exception as e:
                print(e, "Your yt version might have not been patched by me! Retrying with p_size=100.")
                s.annotate_particles(width=(width,width_units), marker = marker, p_size = 100)
        
                    
        if stream_v          : 
            s.annotate_streamlines(("gas","velocity_"+axses[0]),("gas","velocity_"+axses[1]),
                                   density=2,
                                   color="black",
                                   arrowsize=0.,
                                   broken_streamlines=False)
        if stream_m          : 
            s.annotate_streamlines(("gas","magnetic_field_"+axses[0]),("gas","magnetic_field_"+axses[1]),
                                   density=2,
                                   color="midnightblue",
                                   arrowsize=0.,
                                   broken_streamlines=False)
        if hide_axes         : s.hide_axes(draw_frame = True)
        if scale             : s.annotate_scale()
        if grid              : s.annotate_cell_edges()
        if cmap  is not None : s.set_cmap("all", cmap)
        if units is not None : s.set_unit(quant, units)
        if vval  is not None and proj==False : s.set_zlim(quant, zmin = vval.min(), zmax = vval.max())

        return s 

    
        
    def Splotax(self, quant, typefield = "gas", fig = None, kind = "slice", **kwargs):
        """
            This funtion produce a plot using Splot, but it associate it with a matplotlib axis.
            quant : str, quantity to plot
            typefield   : str  , type of the fiel of the quantity
            fig   : figure, a figure to put the plot in 
            kind  : what kind of plot to do. can be "slice" or "proj". Defaults to "slice".
            kwargs: here can be stored arguents to pass to Splot or Pplot.
            return: ax, fig 
        """

        if fig is None: fig = plt.figure(figsize=(8,8))
        
        # initialise the axes
        ax = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                        nrows_ncols = (1, 1),
                        axes_pad = 0.05,
                        label_mode = "L",
                        share_all = True,
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_size="3%",
                        cbar_pad="0%")
        plotter = self.Splot if kind=="slice" else self.Pplot
        s  = plotter(quant, typefield=typefield,**kwargs)#yt.SlicePlot(self.ds,gg,fig=fig, width=width, **kwargs)#, width=(100, 'parsec'))
        gg = (typefield, quant)
        plot        = s.plots[gg]
        plot.figure = fig
        plot.axes   = ax[0].axes
        plot.cax    = ax.cbar_axes[0]
        s._setup_plots()

        return  plot.axes, fig
    
    ###follow particles
    def trace_particles(self, quant, typequant="gas" ):
        """
            This function gives the value of the quantity quant at the location of the particles.
            return: np.array with the requested values
        """
        # sample the field 
        self.ds.add_mesh_sampling_particle_field((typequant, quant), ptype="all")        
        # shouldn't that be sp?
        ad = self.ds.all_data()
        
        ad["all", "cell_index" ]  # Trigger the build of the index of the cell containing the particles
        p_gas = ad["all", "cell_gas_"+quant]
        return p_gas
    
    def particle_data(self, quant = None, only_data = False, cunits = None,  qunits = None):
        """
            # only tested for RAMSES.
            This function compute quantities associated to particles.
            quant    : str ,  particle quantity to extract
            only_data: bool,  extract only data, no coordinates
            cunits   : str ,  units to convert coordinates to 
            qunits   : str ,  units to convert quantity to
        """
        
        if self.sp is None: self.init_sp()
        # extract indexes
        idx = self.sp[('nbody',"particle_index")]
        # sort by indexes 
        p   = np.argsort(idx)
        # if not only_data, extract coordinates
        if not only_data:
            self.xp = self.sp[('nbody', 'particle_position_x')]
            self.yp = self.sp[('nbody', 'particle_position_y')]
            self.zp = self.sp[('nbody', 'particle_position_z')]
            # sort
            self.xp = self.xp[p]
            self.yp = self.yp[p]
            self.zp = self.zp[p]
            # if quant have been defined, extract also quant
            if quant is not None:

                q = self.sp[('nbody', quant)]
                # convert to q units
                if qunits is not None: q=q.to(qunits)
                # sort
                q = np.array(q)
                q = q[p]
                result = [self.xp, self.yp, self.zp, q] if cunits==None else [np.array(self.xp.to(cunits)),
                                                                              np.array(self.yp.to(cunits)),
                                                                              np.array(self.zp.to(cunits)),
                                                                              q]
            else:
                result = [self.xp, self.yp, self.zp] if cunits==None else [ np.array(self.xp.to(cunits)),
                                                                            np.array(self.yp.to(cunits)),
                                                                            np.array(self.zp.to(cunits))]
        else:
            if quant is not None:
                q      =self.sp[('nbody', quant)]
                result =  np.array(q) if qunits==None else np.array(q.to(qunits))
                # sort
                result = result[p]
            else:
                raise Exception("No particle quantity specified. ")

        return result
    

    def multi_part_Splot(self, quant, quantype="gas",
                        fig  = None, width = None,
                        npart           = None,
                        text            = None,
                        columns_rows    = None,
                        text_args       = {"size":10, "color":"black"},
                        tipo            = "square", 
                        orientation     = "vertical",
                        star_index      = True,
                        scale = True, hide_axes = True,
                        vval            = None,
                        scale_text_size = 10,
                        no_cbar         = True,
                        p_size          = ("mass", "Msun",100),
                        marker          = r"$\odot$",
                        **kwargs):
        """
            This function produces multiple slice plots centered on a number of particles.
                    
            This function produces Slices of a simulation output.
            quant       : str,            quantity to plot
            fig         : figure,         a figure to put the plot in 
            width       : float,          width of the region to plot
            npart       : int|range|list, gives indexes or number of particles to print
            text        : str,             what to write on the plots
            columns_row : tuple(int,int), how many columns and rows
            text_args   : dict,           arguments to pass to the call to plt.text in yt.sliceplot             
            tipo        : str,            how to auto check for columns and rows
            orientation : str,            in case of autocheck, set what will be the long side
            star_index  : bool,           print or not the indexes of the star
            scale       : bool,           annotate a scale bar in the last plot
            hide_axes   : hide aces
            vval        : tuple|list|np.array, select only a range of values.
            scale_text_size: int,        fontsize of the text for the scale bar
            no_cbar     : bool,          hide color bars
            p_size      : float|int|(field, units, int|float|func): give rules to scale the size of the particle markers
                        THIS ACTUALLY REQUIRES TO PATH THE YT FUNTION
            marker      : str,            select the marker for the particles positions 
            
            kwargs      : here can be stored arguents to pass to Splot.
        """
        
        if fig is  None  : fig   = plt.figure(figsize = (8,8))
        if width is  None: width = float(self.width[0].d)

        
        gg        = quant
        # save positions 
        positions = self.particle_data()
        # save mass of the particles
        mass      = self.particle_data(quant = "particle_mass", only_data = True, qunits = "Msun")
        
        # define how many and/or wich particle make the plots about 
        if npart is None:
            # case all particles
            npart   = len(positions[0])    
            loopvar = np.array(range(0, npart))
        elif isinstance(npart, list) or isinstance(npart, range):
            # selected particles
            loopvar = np.array(npart)
            npart   = len(loopvar)
        elif isinstance(npart,int):
            # just up to a number of particle  
            loopvar = np.array([npart])
            npart   = len(loopvar)
        else: 
            raise("npart type not acceptable.")
        
        # set how many rows
        if columns_rows is None:

            # this function is in other
            columns_rows = compute_rows_columns(len(loopvar), tipo, orientation)
        if vval is False and len(loopvar) > 1 and not no_cbar:
            # in this case, put a colorbar for each plot
            cbar_mode = "each"  
        elif vval is False and len(loopvar)>1 and no_cbar:
            # no colorbars at all
            cbar_mode = None  
        else:
            # just one colorbar for all of them 
            cbar_mode = "single"
        # define axes 
        ax = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                        nrows_ncols  = columns_rows,#(3, 5),
                        axes_pad     = 0.05,
                        label_mode   = "L",
                        share_all    = True,
                        cbar_location= "right",
                        cbar_mode    = cbar_mode,
                        cbar_size    = "3%",
                        cbar_pad     = "0%"
                        )
        if vval is None:
            # find auto colorbar limits  
            if self.sp is None: self.init_sp()
            vval = self.sp[(quantype, quant)]
            vval = np.array([vval.min(),vval.max()])
        elif vval is False:
            # ignore limits and set vval to false
            vval = None
        else: 
            # take vval as supplied
            vval = np.array(vval)
            vval = np.array([vval.min(),vval.max()])
        
        nonetext_check = False
        for i, j in zip(loopvar, range(0,len(loopvar))):  
            if text is None or nonetext_check: 
                # auto 
                mass  = self.particle_data(quant="particle_mass",only_data=True, qunits="Msun")#[i]
                texpl = self.particle_data(quant="particle_birth_time",only_data=True, qunits="Myr")#[i]
                text1 = r" M = $%.2f \,\, \rm M_\odot$ "%(mass[i])
                text2 = r"$\tau_{SN} = %.2f\,\, \rm Myr$ "%(texpl[i])
                text  = text1+"\n"+text2
                nonetext_check=True
            elif text is False:
                # none
                text  = ""
            else: 
                # supplied
                pass
            # find centers of each plot
            centrei = self.findpos(sp=self.sp ,poss=positions, part=i)
            # produce slice 
            s = self.Splot(gg, width = width, particles = True, center = centrei, time = False, hide_axes = hide_axes,
                          vval   = vval,
                          p_size =p_size,
                          marker =marker,
                         **kwargs)
            
            # annotate scale in the last plot only 
            if j == len(loopvar)-1: s.annotate_scale(text_args={'size':scale_text_size})#text_args=text_args)
            # annotate text
            s.annotate_text((0.1, 0.8), text, coord_system="axis", text_args=text_args)
            # annotate star indexes
            if star_index: s.annotate_text((0.1, 0.1), r"$\bigstar = %d$"%i, coord_system="axis", text_args=text_args)
            
            plot        = s.plots[gg]
            plot.figure = fig
            plot.axes   = ax[j].axes
            plot.cax    = ax.cbar_axes[0 if cbar_mode=="single" else j ]
            
            s._setup_plots()

        return  plot.axes, fig
#### 3D PLOTS
    def plot3D(self,
                quant           ,#= "magnetic_field_magnitude",
                outpath         = None,
                threshold_field = "velocity_magnitude",
                axes            = ["x","y","z"],
                load_data       = True,
                save_data       = True,
                threshold       = 0., 
                scale_axes      = [costs.pc,costs.pc,costs.pc],
                scaleq          = None,
                operation       = ">",
                clim            = None,
                pbounds         = None, 
                point_size      = 0.8,
                split=None,
                render_points_as_spheres = False, 
                log = True, 
                log_axes = [False, False, False], 
                #clip_axes=
                log_func = np.log10, 
                diffuse  = None, 
                tf       = None, 
                cmap     = "rainbow",
                particles= True,
                trajectories_from_file = True,
                trajectories_path      = "./trajectories",
                plotter  = None,
                labels   = None, 
                density_phase= False,
                centered = True,
                gaussian=False, 
                rescale=-1
                ):
        if tf is None: tf = tff.tff(model="sin", alpha=0)
        if not isinstance(threshold_field,list): threshold_field = [threshold_field]
        if not isinstance(threshold, list): threshold=[threshold]
        if not isinstance(operation, list): operation=[operation]
        if plotter is None: plotter = Plotter(off_screen = False)
        if outpath is None: outpath = self.directory
        
        
        
        if isinstance(tf,tuple):
            tf = tff.tff( model = tf[0], i = tf[1])
        elif isinstance(tf,str):
            if tf=="linear":
                pass
            else:
                tf = tff.tff( model = tf)
        elif isinstance(tf, float) or isinstance(tf, int) :
            tf = tff.tff( model = tf)
        elif isinstance(tf, np.ndarray):
            """ just keep it as it is.  """
            pass
        
        else:
            ## so if it is a function basically 
            tf = tff.tff( model = tf)
        
        if particles:
            if "x" in axes and "y" in axes and "z" in axes:
                
                x, y ,z = self.particle_data(cunits="pc")
                if centered: 
                    x-=0.5*(self.width.to("pc").d)[0]
                    y-=0.5*(self.width.to("pc").d)[0]
                    z-=0.5*(self.width.to("pc").d)[0]
                points = np.column_stack((x, y, z))
                point_cloud = PolyData(points)

                plotter.add_points(
                            point_cloud,
                            #clim=clim, 
                            #scalars=qname,
                            point_size = 10,#0.8,#0.8,
                            cmap     = cmap,
                            style    = "points_gaussian",
                            emissive=True,
                            render_points_as_spheres = render_points_as_spheres,#render_points_as_spheres,
                            diffuse=1
                            )
                if trajectories_from_file:
                    self.plot_trajectories_from_files(plotter=plotter,negative_time=True, data_path=trajectories_path,centered=centered)
            else:
                print("Trajectories not allowed with these axes: %r"%axes)
        if labels is None: labels=axes.copy()
        qname = quant 
        mags, coords = utils_3D.load_and_save_simple(self,
                                            quant           = quant,
                                            rs_path         = outpath,
                                            threshold_field = threshold_field,
                                            axes            = axes,
                                            load_data       = load_data,
                                            save_data       = save_data,
                                            )
        threshold_field, quant =  mags
        #print(mags, coords)
        #print(mags)
        
        
         
        mags, coords = utils_3D.setup_3Ddata(
                mags,
                coords,
                split=split,
                scale_axes= scale_axes,
                scaleq    = scaleq,
                threshold = threshold,
                operation = operation,
                clim      = clim,
                centered=centered, 
                box=(self.width.to("pc").d)[0]
                )
        threshold_field, quant =  mags
        del threshold_field   
        
        if log_axes is not None: 
            for i,log_scale in enumerate(log_axes):
                #print(coords[i])
                if log_scale: coords[i] = log_func(coords[i])
                
        if density_phase:
            
            h,edges = utils_3D.bins3d(coords[0],
                                    coords[1],
                                    coords[2],
                                    bins=512,
                                    density=False
                                    ) 
            qname = "phase distr"
            quant = h.flatten()
            
            quant /= quant.max()
            
            x,y,z = np.meshgrid(edges[0],edges[1],edges[2])
            x,y,z =[x.flatten(),y.flatten(),z.flatten()]
            quant, coords =utils_3D.clip_arrays([x,y,z], quant, 0., ">")
            
        #print(len(quant), len(coords[0]),len(coords[1]), len(coords[2]))
        
        if diffuse=="adaptive": diffuse = utils_3D.set_diffusive(quant)
        if point_size=="adaptive": point_size=utils_3D.set_size(quant, 3, 1.5)
        #print(diffuse)
        #sys.exit()
        plotter = utils_3D.plot3D(coords, 
                            qname, 
                            quant,
                            gaussian=gaussian, 
                            pbounds    = pbounds, 
                            point_size = point_size,
                            render_points_as_spheres = render_points_as_spheres, 
                            log        = log, 
                            diffuse    = diffuse, 
                            tf         = tf, 
                            cmap       = cmap,
                            plotter    = plotter,
                            labels     = labels
                            )
        
  

        
        
        return plotter
    #@staticmethod        
    def plot_trajectories_from_files(
                        self,
                        scale_l   = costs.pc,
                        scale_t   = costs.myr,
                        negative_time= False,
                        #dense     = False,
                        plotter   = None,
                        data_path = None,
                        base_names= "snp",
                        fmt       = "%02d",
                        ext       = ".dat",
                        centered=True):
        from pyvista import    MultipleLines
        if plotter is None: 
            
            plotter = Plotter()
        
        if data_path is None: data_path = self.directory +"/../"+ "trajectories/"
        data_path +="/"  
        x, y ,z, indx = self.particle_data(quant="particle_index", cunits="pc")
        #masses = self.particle_data(only_data=True, quant="particle_mass", qunits="Msun")

        #first = True
        current_t = self.time.to("s").d/scale_t
        #time = np.linspace(0,self.time.to("s").d/scale_t, 1000,endpoint=True)
        for idx in indx:
            tpath = data_path +  base_names + fmt%(int(idx)-1) + ext 
            if not os.path.exists(tpath): raise Exception("No trajectories found.")
            with open(file=tpath, mode="r") as f:
                data = np.loadtxt(f)
            t = data[:,0] * costs.myr / scale_t
            
            if negative_time: t = t - t[0] 
            x = data[:,1] * costs.pc / scale_l
            y = data[:,2] * costs.pc / scale_l
            z = data[:,3] * costs.pc / scale_l
            if not centered:
                x +=  self.width[0].d*0.5
                y +=  self.width[0].d*0.5
                z +=  self.width[0].d*0.5
            points_set = [[xi,yi,zi] for (xi,yi,zi) in zip(x[t<= current_t],
                                                    y[t<= current_t],
                                                    z[t<= current_t]) ]
            mesh = MultipleLines(points=points_set)
            
                
            plotter.add_mesh(mesh, color='w', line_width=1)
        #plotter.show()
        
        return plotter        
#ds.Splotax(quant)


    def streamlines3D(self,
                    quant           = None,
                    level           = 7,
                    log             = True,
                    outpath         = None,
                    key             = "magnetic_field",
                    bounds          = None,
                    pbounds         = None,
                    threshold_field = "velocity_magnitude",
                    threshold       = 5e5,
                    num_points      = 500,
                    from_source     = True, 
                    #stream_opacity= None,
                    tf              =None,  
                    max_integration_time  =1000,
                    cut_below       = True,
                    cmap            = "rainbow", 
                    clim            = None,
                    selected_points = None,
                    clip            = False,
                    load_grid       = False,
                    load_data       = True,
                    save_data       = True,
                    save_grid       = True,
                    plotter         = None,
                    centered= True
                    ):
        if plotter is None: plotter = Plotter(off_screen = False)
        #print(plotter)
        #sys.exit()
        if tf is None: tf= tff.tff(model="sin", alpha=0)
        if outpath is None: outpath = self.directory
        if isinstance(tf,tuple):
            tf = tff.tff( model = tf[0], i = tf[1])
        elif isinstance(tf,str): 
            tf = tff.tff( model = tf)
        elif isinstance(tf, float) or isinstance(tf, int) :
            tf = tff.tff( model = tf)
        elif isinstance(tf, np.ndarray):
            """ just keep it as it is.  """
            pass
        else:
            ## so if it is a function basically 
            tf = tff.tff( model = tf)
            
            
            
        if quant is None: quant = key+"_magnitude" 
        LOG = log
        ds  = self      
        thresh_mag_data, vec_data = utils_3D.load_and_save_data(
                                            ds             = ds,
                                            level          = level,
                                            rs_path        = outpath,
                                            quant          = quant,
                                            load_data      = load_data,
                                            save_data      = save_data,
                                            threshold_field = threshold_field,
                                            vectorial_data = True,
                                            key            = key,
                                            )



        #vec_data = magnitude_cut, thresh_data_cut, vec_x_cut, vec_y_cut, vec_z_cut  
        grid,thresh_mag_data, vec_data = utils_3D.extract_cut(
                                            ds.width[0],
                                            thresh_mag_data,
                                            vec_data,
                                            bounds = bounds)

        gridcut, quant_name = utils_3D.load_and_save_grid(
                                            level=level, 
                                            rs_path         = outpath,
                                            grid            = grid,
                                            key             = key,#quant_name="mmm"
                                            thresh_mag_data = thresh_mag_data,
                                            LOG             = LOG,
                                            vec_data  = vec_data,
                                            load_grid = load_grid,
                                            save_grid = save_grid)
        
        if centered: 
            cnt=tuple([-self.width[0].to("pc")/2]*3)
            gridcut=gridcut.translate(cnt, inplace=False)
            #print(gridcut,gridcut.center)
            #sys.exit()
            
        plotter = utils_3D.streamplot(
                            gridcut               = gridcut,
                            thresh_mag_data       = thresh_mag_data,
                            quant_name            = quant_name,
                            key                   = key,
                            LOG                   = LOG,
                            threshold             = threshold,
                            num_points            = num_points,
                            from_source           = from_source,
                            tf                    = tf,
                            max_integration_time  = max_integration_time,
                            cut_below             = cut_below,
                            cmap                  = cmap, 
                            clim                  = clim,
                            selected_points       = selected_points,
                            clip                  = clip,
                            pbounds               = pbounds,
                            plotter               = plotter,
                            )
        
        #plotter.show(auto_close=False)
        return plotter





#### 1D PLOT
    def oneplot(self, 
                x, y, w = None,
                binned  = False,
                bins    = 128,
                maxrange= None,
                fig     = None , ax     = None,
                xtype   = "gas", ytype  = "gas", wtype = "gas",
                logx    = False, logy   = True,
                x0 = 0., y0=0., 
                xlabel  = None , ylabel = None,
                scalex  = 1    , scaley = 1,
                plotstyle = "plot",
                return_data = False, 
                **kwargs,
                ):
        """
            This function produces 1D plots. 
            x              : str,    name of the field on the x axis 
            y              :              "                y axis
            w              : str,         "                to be used as weight,
                            only used when binned = True
            binned         : bool, use the funtion histogram to compute averaged distributions.
            bins           : int, number of bins, only used if binned is true.  
            maxrange       : np.array(x0,x1), max range for the x axis and the binning. Only used when 'binned' is True
            fig            : figure, if no figure is supplied, create one
            xtype / ytype / wtype : str,    types of the field to plot on
            logx  / logy   : bool,   use log scale on the relative axis
            xlable/ylabel  : str,    label for the relative axis 
            scalex/scaley  : str|int|float, scalte to conver the values
            plotstyle      : str, style of the plot. Can be either "plot" or "scatter"
            return_data    : bool, return data too
            return: ax, fig, if return data is true ax, fig, x, y
        """
        if w is None and binned is True: raise Exception("'w' must be specified if 'binned' is True.")
        if self.sp is None: self.init_sp()
        xname , yname = x , y
        x = self.sp[(xtype, x)]
        y = self.sp[(ytype, y)]
        if w is not None: w = self.sp[(wtype, w)]
        if  isinstance(scalex, str):
            x  = x.to(scalex)
            ll = r"%s [$\rm{%s}$]"%(xname, x.units)
        else:
            
            if scalex == 1:
                ll = r"%s [$\rm{%s}$]"%(xname, x.units) 
            else:
                ll = r"%s [$\rm{%s}$/%.3E]"%(xname, x.units, scalex) 
                x /= scalex

        ax.set_xlabel(ll if xlabel is None else xlabel )
        
        if  isinstance(scaley, str):
            y  = y.to(scaley)
            ll = r"%s [$\rm{%s}$]"%(yname, y.units) 
        else:
            
            if scaley == 1 or scaley == 1. :
                ll = r"%s [$\rm{%s}$]"%(yname, y.units) 
            else:
                ll = r"%s [$\rm{%s}$/%.3E]"%(yname, y.units, scaley)
                y /= scaley 
        ax.set_ylabel(ll if ylabel is None else ylabel )        
        p = np.argsort(x)
        # sort the data in the x axis and the y data according to that
        x = x[p]
        y = y[p]
        #if w is not None : w = w[p]
        if binned: 
            x, y = utils.oneD_avg_weighted_bins(x.d,y.d, w.d, bins = bins, maxrange = maxrange)
        x -= x0
        
        y -= y0
        
        if plotstyle == "plot":
           ax.plot(x, y, **kwargs)
        elif plotstyle == "scatter":
           ax.scatter(x,y, **kwargs)
        else:
            
            raise Exception("This plot style is not available. You can only do 'plot' or 'scatter' ")
        
        if logx: ax.set_xscale("log")
        if logy: ax.set_yscale("log")
        if "label" in kwargs: ax.legend()
        fig.tight_layout()

        return [ax , fig, x, y]  if return_data else [ ax , fig]
    
    
    def cached_sp(self, key):
        cachedir = "sp_cache"
        if not os.path.exists(f"{self.directory}/{cachedir}"): os.mkdir(f"{self.directory}/{cachedir}")
        if self.cached:
            try: 
                data = pickle.load(open(f"{self.directory}/{cachedir}/{key}.pkl", "rb"))
                print("Loaded from cache %s"%str(key))
            except:
                print("Not in cache %s, reading."%str(key))
                data = self.sp[key]
                pickle.dump(data, open(f"{self.directory}/{cachedir}/{key}.pkl", "wb"))
        else:
            data = self.sp[key]

        return data

###phase plots
    def phase(self, x, y, w = "volume",
               xtype     = "gas", ytype = "gas", weigthtype = "gas",
               bins      = 512,
               cmap      = "hsv",
               cmaplabel = None,
               xlabel    = None,
               ylabel    = None,
               fig = None, ax = None, 
               scatter=False,
               logx = True, logy = True, logw = True, sumist = False,
               absx=False,absy=False,
               only_scalar=False,
               subregions = None,
               contours  = False,
               **kwargs):
        """
            This function prints phase plots of the given quantities.
            x,y,w                   : str, quantities to be plotted in the x and y anxis and to be used as weights.
            xtype, ytype, weightype : str, types of the fields
            absx,absy               : bools, use absolute values.
               
            bins                    : int, how many bins
            cmap                    : str, what cmap shoul be used
            x/ylabel/cmaplabel      : str, labels of the various axes and colorbar
            fig, ax                 : existing fig and axes that can be used to plot into 
            logx/y/w                : bool, use log scale for the given quantity
            return ax, fig
        """
        if isinstance(contours, str): countours = eval(contours)

        # autoset labes
        if xlabel    is None: xlabel = fr"{str(myunits.symbols[x])}"+ str(myunits.units[x]) #paper.units.labels[x]
        if ylabel    is None: ylabel = fr"{str(myunits.symbols[y])}"+ str(myunits.units[y])
        if cmaplabel is None:
            cmaplabel = r"$\mathrm{Pdf}$" if w is None else  fr"{str(myunits.symbols[w])}"+ str(myunits.units[w])
        else:
            cmaplabel = r"%s"%cmaplabel
        if self.sp   is None: self.init_sp()


        # generate fig and ax if not supploes
        fig = plt.figure() if fig is None else fig
        ax  = fig.add_subplot(111) if ax is None else ax
        
        # extract values 
        # x = self.sp[(xtype,x)].in_cgs().value
        # y = self.sp[(ytype,y)].in_cgs().value
        xname = str(x)
        x = self.cached_sp((xtype,x))
        yname = str(y)
        y = self.cached_sp((ytype,y))
        if yname=="temperature": y*= costs.mu_mol


        if w is not None:
            wname = str(w)
            # if wname=="mass":
                # w = self.sp[(weigthtype, w)].to("Msun").value 
            # else:
            w = self.cached_sp((weigthtype, w)) 
        else:
            w = None
        if only_scalar:
            scalar = self.cached_sp(("ramses", "hydro_scalar_00"))
            scalar2 = self.cached_sp(("gas", "velocity_magnitude"))


            tresh, tresh2 = 1e-10, 0#15e5
            
            cond = np.logical_and(scalar>tresh,scalar2 > tresh2)
        else: 
            cond = np.ones_like(x, dtype=bool)
        if absx:
            x = abs(x)
        if absy:
            y = abs(y)
        if logx:
            x = np.log10(x)
        if logy:
            y = np.log10(y)
        
        if logw and w is not None: 
            w = np.log10(w)
        # compute histograms  and make plots
        if scatter:
            im=ax.scatter(10**x,10**y,
                        c = 10**w,
                        cmap=cmap,
                        **kwargs,
                        norm=colors.LogNorm() if logw else None)
            if logx: ax.set_xscale("log")
            if logy: ax.set_yscale("log")
            fig.colorbar(im, ax = ax)
        else:
            N = 0
            ggs = {0:"cond"}#{0:"not", 1:"cond"}
            
            for gg in ggs.values():
                if gg == "cond":
                    condt = cond
                    
                else:
                    condt = np.logical_not(cond)

                xO = x[condt]
                yO = y[condt]
                wO = w[condt] if w is not None else None

                # if gg == "cond":
                #     cmap = plt.cm.get_cmap("Greys")
                # else:
                #     cmap = plt.cm.get_cmap("Greens")

                    
                if wO is not None:
                   
                    sums, xbins, ybins = np.histogram2d(xO, yO, bins = bins)#, range=[[-26, -21.5], [2, 7.5]])
                    counts,xbins, ybins= np.histogram2d(xO, yO, bins = [xbins, ybins],
                                                        weights =10**wO if logw else wO)
                    # print(xbins, ybins)
                    xbins = 0.5 * (xbins[1:] + xbins[:-1])
                    ybins = 0.5 * (ybins[1:] + ybins[:-1])
                    if logx: xbins   = 10**xbins
                    if logy: ybins   = 10**ybins
                    counts = (counts.T)
                    
                    with np.errstate(divide='ignore', invalid='ignore'):
                        if not sumist: 
                            sums   = (sums.T  )
                            goin = counts/sums
                            
                            img = ax.pcolormesh(xbins, ybins, goin , cmap = cmap,
                                                norm=colors.LogNorm() if logw else "linear",
                                                rasterized=True,
                                                **kwargs)
                        else:
                            print("-- sum, does not work good. ")
                            ccond  = counts < 1.0
                            counts = np.where(ccond, np.nan, counts)
                            if wname == "mass": 
                                goin = counts / 1.989e33
                            datatodump = [goin, xbins, ybins]
                            pickle.dump(datatodump, open(f"histdump_{xname}_{yname}_{wname}_{self.pathinfo["outnumb"]:05d}.pk", "wb"))
                            img = ax.pcolormesh(xbins, ybins, goin, 
                                                norm=colors.LogNorm() if logw else "linear",
                                                cmap = cmap,
                                                rasterized=True,
                                                
                                                **kwargs)
                            if gg == "not": continue
                            
                            goin[np.isnan(goin)] = 0
                            # Calculate the total mass in the entire phase space
                            total_mass_all = np.nansum(goin)  # Sum all non-NaN values

                            # Place the total mass text in the top right corner
                            if only_scalar:
                            
                                Y0TEXT = 0.15#0.95
                                TEXT = r"$\mathrm{M}_{\rm SNR} = %.2f\,\,\rm{M_\odot}$"
                            else:
                                Y0TEXT = 0.95
                                TEXT = r"$M_{\rm tot} = %.2f\,\,\rm{M_\odot}$"
                            
                            ax.text(0.95, Y0TEXT, TEXT % total_mass_all, 
                                ha='right', va='top', transform=ax.transAxes, 
                                color='black', fontsize=8, bbox=None)
                            


                            loggoin = np.log10(goin)
                            

                            mylevels = [-2, -1.7]
                            mylevels = [-2, -1.5]
                            if contours:  
                                cs = ax.contour(xbins, ybins, gaussian_filter(loggoin,1), levels = mylevels,
                                                                #alpha=0.1, filled = True
                                                                linestyles=["dashed","solid"],
                                                                colors ="black",#"red",
                                                                linewidths = 0.2,
                                                                )  # alpha=0 to hide fill if desired
                               
                                #ax.clabel(cs, cs.levels, fmt=fmt, fontsize=4)
                            subregions = None
                            # subregions =  [  [[1e2,2e3],[2e-25,5e-23]],
                            #                  [[2e3,1e4],[2e-25,1.4e-24]],
                            #                  [[2e3,3.5e4],[9.5e-25,2e-22]],
                            #                  [[3.5e4,3e5],[2e-25,1e-23]],
                            #                  [[3e5,2e7],[1e-26,3e-24]]]
                            colorstmp = ['red', 'blue', 'green', 'purple', 'orange', 'cyan', 'magenta', 'yellow']
                            
                            
                            if subregions is not None:

                                    
                                yposs = [Y0TEXT-0.05*(i+1) for i in range(len(subregions))] 
                                        
                                X, Y = np.meshgrid(xbins, ybins)
                                
                                # Define a small multiplicative factor for the offset
                                offset_factor = 0.02  # Adjust this value as needed (e.g., 0.05 for a 5% offset)
                                total_mass = 0.0
                                excluded_mask = np.zeros_like(goin, dtype=bool)

                                for isub, sub in enumerate(subregions):
                                    
                                    Tlims, rholims = sub
                                    T_min, T_max = Tlims
                                    rho_min, rho_max = rholims
                                    if False:
                                        T_min_adj = T_min * (1 + offset_factor)
                                        T_max_adj = T_max * (1 - offset_factor)
                                        rho_min_adj = rho_min * (1 + offset_factor*0.75)
                                        rho_max_adj = rho_max * (1 - offset_factor*0.75)
                                        rect = Rectangle(
                                            (rho_min_adj, T_min_adj),  # Adjusted lower-left corner
                                            rho_max_adj - rho_min_adj,  # Adjusted width
                                            T_max_adj - T_min_adj,      # Adjusted height
                                            fill=False,                 # No fill, just an outline
                                            edgecolor=colorstmp[isub],
                                            linewidth=0.5,
                                            alpha=0.5
                                        )
                                        ax.add_patch(rect)
                                    # Use np.logical_and to get proper boolean arrays for the current subregion
                                    cond_T = np.logical_and(Y >= T_min, Y <= T_max)
                                    cond_rho = np.logical_and(X >= rho_min, X <= rho_max)
                                    current_mask = np.logical_and(cond_T, cond_rho)

                                    # Restrict the excluded_mask to the current subregion's bounds
                                    subregion_mask = np.zeros_like(goin, dtype=bool)
                                    subregion_mask[np.ix_(
                                        np.logical_and(ybins >= T_min, ybins <= T_max),
                                        np.logical_and(xbins >= rho_min, xbins <= rho_max)
                                    )] = True

                                    # Exclude data already included in previous subregions
                                    current_mask = np.logical_and(current_mask, ~excluded_mask)

                                    # Update the cumulative excluded mask
                                    excluded_mask = np.logical_or(excluded_mask, current_mask)

                                    # Get the 1D indices for xbins and ybins that satisfy the condition
                                    idx_rho = np.where((xbins >= rho_min) & (xbins <= rho_max))[0]
                                    idx_T = np.where((ybins >= T_min) & (ybins <= T_max))[0]

                                    if len(idx_rho) == 0 or len(idx_T) == 0:
                                        print("No data in this region.")
                                        continue

                                    # Extract the subarray from goin using np.ix_ for proper 2D indexing
                                    xbinst = xbins[idx_rho]
                                    ybinst = ybins[idx_T]
                                    goint = goin[np.ix_(idx_T, idx_rho)]

                                    # Apply the current mask to the subarray
                                    goint = np.where(current_mask[np.ix_(idx_T, idx_rho)], goint, np.nan)
                                    # Logarithmic transformation for contour plotting
                                    loggoint = np.log10(goint)#, where=(goint > 0))
                                    
                                    
                                    ax.contour(xbinst, ybinst, loggoint, levels = mylevels, 
                                                                alpha=1, colors=colorstmp[isub], linestyles="solid", linewidths=0.2)  # alpha=0 to hide fill if desired

                                    # Create meshgrid for coordinate mapping
                                    
                                    gointosumt = goint.copy()

                                    # Iterate through each filled contour region
                                    xpos = 0.95
                                    ypos = yposs[isub]

                                    tmpmass = np.nansum(gointosumt)
                                    total_mass += tmpmass
                                    ax.text(xpos, ypos, r"$%.2f\,\,\rm{M_\odot}$"%tmpmass,
                                            color=colorstmp[isub] , fontsize=8, ha='right', va='top', transform=ax.transAxes,
                                            bbox=None,
                                            
                                            )
                                    print(f"In subregion {isub}, {sub}, mass is {tmpmass:.2f}, of a total {total_mass:.2f}")

                                    N+=1
                                    
                else:
                    
                    sums, xbins, ybins   = np.histogram2d(xO, yO, bins=bins)
                    counts, xbins, ybins = np.histogram2d(xO, yO, bins=[xbins, ybins])#, density = True)
                    xbins = 0.5 * (xbins[1:] + xbins[:-1])
                    ybins = 0.5 * (ybins[1:] + ybins[:-1])
                    
                    
                    if logx: xbins   = 10**xbins
                    if logy: ybins   = 10**ybins
                    
                    
                    # suppress possible divide-by-zero warnings
                    with np.errstate(divide='ignore', invalid='ignore'):
                        if logw:
                            img = ax.pcolormesh(xbins, ybins, counts.T  ,
                                                cmap = cmap, norm=colors.LogNorm(),**kwargs) #**kwargs)
                        else: 
                            # only color pixels with at least 1 cell count
                            counts=counts.T
                            ccond = counts < 1.0
                            counts= np.where(ccond,np.nan, counts)
                            
            
                            img = ax.pcolormesh(xbins, ybins, counts  ,
                                                cmap = cmap, **kwargs)
                #new_ax =  inset_axes(ax, width="70%", height="5%", loc="upper center") 
            divider = make_axes_locatable(ax)
            cax = divider.new_vertical(size = '5%', pad=0.1, )
            fig.add_axes(cax)
            cb     = fig.colorbar(img,
                                    #ax = cax,
                                    cax=cax,
                                    extend = 'both',
                                    orientation = 'horizontal',
                                    location = 'top'
                                    )
            cb.set_label(label = cmaplabel, **paper.figfontstyle)
            
            if logx: ax.set_xscale("log")
            if logy: ax.set_yscale("log")
            ## set scalar formatter if there isn't enough range
            if abs(x.max()-x.min()) < 1:
                ax.xaxis.set_major_formatter(ScalarFormatter())
                ax.xaxis.set_minor_formatter(ScalarFormatter())
            if abs(y.max()-y.min()) < 1:
                ax.yaxis.set_major_formatter(ScalarFormatter())
                ax.yaxis.set_minor_formatter(ScalarFormatter())
            		
   #cb.ax.xaxis.set_ticks_position('top')
            #cb.ax.xaxis.set_label_position('top')
			
            #cbar = fig.colorbar(img, ax = ax)
            
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)

        return ax, fig

###check turbulence
    
    def spectrum(self, level, workers = None, key = "vel", filter = None, 
                     ax = None, fig = None, plot = True, color = "black", label = "", kol = True, lw = 5, lwp = 2.5, interface = scipy,
                     colorperp = "red", colorpara="midnightblue",
                     labelperp="Perp.", labelpara="Para." ):
        
        # if plot is needed 
        if plot:
            fig = plt.figure() if fig is None else fig
            ax  = fig.add_subplot(1,1,1) if ax is None else ax
        # only quantities to be analised are magnetic field and velociy
        if key not in ["mag","vel" ]: sys.exit("The key '%s'  is not included in the available ones."%key) 
        
        # selecting fields to exctract data of
        if key == "vel": base="velocity_"
        if key == "mag": base="magnetic_field_"
        axis = "xyz"
        fieldlist = [("gas", base+dir) for dir in axis ] 
        #fieldlist=[]
        #for dir in axis: fieldlist += [("gas", base+dir)] 

        # building regular grid
        
        low       = self.ds.domain_left_edge
        dims      = self.ds.domain_dimensions * 2**level
        

        # extract the fields data with smoothed grid 

        cube = self.ds.smoothed_covering_grid(level, left_edge = low, dims = dims, fields = fieldlist)
        # what is this...?

        f = cube[filter[1]].d if filter is not None else None
        
        nx, ny, nz = dims
        Kk = np.zeros((nx // 2 + 1, ny // 2 + 1, nz // 2 + 1))

        # compute fast fourier transform  for each component of the vector field and sum over Kk
        for vel in fieldlist:
            Kk += 0.5 * fft_comp(self.ds,
                                      cube,
                                      vel,
                                      interface,
                                      workers, 
                                      filter = filter,
                                      f=f)

        # wavenumbers

        L = (self.ds.domain_right_edge - self.ds.domain_left_edge).d
        kx = interface.fft.rfftfreq(nx) * nx / L[0]
        ky = interface.fft.rfftfreq(ny) * ny / L[1]
        kz = interface.fft.rfftfreq(nz) * nz / L[2]
        
        # physical limits to the wavenumbers
        kmin = np.min(1.0 / L)
        kmax = np.min(0.5 * dims / L)
        # binning the wavenumbers
        kbins = np.arange(kmin, kmax, kmin)
        N = len(kbins)

        # bin the Fourier KE into radial kbins

        kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
        
        k2 = kx3d**2 + ky3d**2 + kz3d**2  
        k = np.sqrt(k2)   
        if key == "mag":
            # direction of the mean magnetic field, if any. 
            reference_direction  = np.array([0,1,0])
            k2r = kx3d**2*reference_direction[0] + ky3d**2*reference_direction[1] + kz3d**2*reference_direction[2]
            
            # compute sum over the parallel wavenumbers
            k_parallel = np.sqrt(k2r)     

            # compute the perpenficular one
            k_perpendicular = np.sqrt(k**2 - k_parallel**2)
            
            # paraLlel 
            kpara, E_para,_,_ = calc_spectrum_E(k,Kk, k_parallel, kbins, N)
            
            ##perp

            kperp, E_perp,_,_ = calc_spectrum_E(k,Kk, k_perpendicular, kbins, N)

        k, E_spectrum, kmax, Emax = calc_spectrum_E(k,Kk, k, kbins, N)
        
        if plot:
            if  key=="mag":
                ax.loglog(kperp, E_perp, label=labelperp,lw=lwp , color=colorperp)
                ax.loglog(kpara, E_para, label=labelpara,lw=lwp , color=colorpara)
                #ax.loglog(k, E_para+E_perp, label="tot")
                

            ax.loglog(k, E_spectrum, lw = lw, color = color, label = label)
            if kol: ax.loglog(k, Emax * (k / kmax) ** (-5.0 / 3.0), ls=":", color="0.5")
            #ax.set_xscale("log")
                #ax.set_yscale("log")
            ax.set_xlabel(r"k")#r"$\lambda(k)\,\,[pc]$")
            ax.set_ylabel(r"$E(k)dk$ ")
            ax.legend()
        del cube, f
        return  [k, E_spectrum, ax, fig] if plot else [k, E_spectrum]
    

    
    def hist(self, quant, 
        w = None,
        fig = None, ax = None,
        typequant = "gas", 
        typew     = "gas", 
        scaleq    = 1., scalew = 1.,
        invw      = False,
        bins      = 50,
        #density   = False,
        #label     = "",
        stat      = True,
        scolor     = "midnightblue",
        salpha     = 1.,
        absolute  = False,
        absolutew = False,
        #log = True,
        logx = True,
        xlabel = None, ylabel = None,
        **kwargs
        ):
        """ 
            This funtion procudes histograms of a given quantity from a simulation snapshots.
            quant       : str  , field to plot the histogram of 
            w           : str  , field to be used as a weight 
            typequant   : str  , type of the field of quant
            typew       : str  , type of the field of w
            scaleq      : float, scaling factor for q 
            scalew      : float,        "           w
            invw        : bool ,
            bins        : int  , number of bins to use in the histogram
            stat        : bool , compute and plot mean and std of quant 
            scolor      : str  , color of the stat plot
            salpha      : float, alpha of the stat plot
            absolute    : bool , use absolute values for quant
            absolutew   : bool ,            "            w
            logx        : bool , set the x axis as a log
        xlabel, ylabel  : str(s), labels
            return: [ax, mean, std] if stat else ax
        """
        # set the labels
        if xlabel is None: xlabel = quant       
        if ylabel is None:  
            # check is density is in the supplied kwargs
            if 'density' in kwargs:
                density  = kwargs["density"]
            else:
                density = False # this is default in plt.hist()
            ylabel = ("N" if not density else "p") if w == None else w

        # init dataset if needed 
        if self.sp is None: self.init_sp()
    
        # create figure and axis if they are not supplied
        fig = plt.figure() if fig is  None else fig
        ax  = fig.add_subplot(111) if ax is None else ax

        # extract values to make histogram 
        if isinstance(scaleq, float):
            quant = (self.sp[(typequant, quant)].in_cgs()).v/scaleq
        elif isinstance(scaleq, str):
            quant = (self.sp[(typequant, quant)].to(scaleq)).v

        # set weight
        if not ( w is None): 
            if isinstance(scaleq, float):
                w = (self.sp[(typew, w)].in_cgs()).v/scalew
            elif isinstance(scaleq, str):
                w = (self.sp[(typew, w)].to(scalew)).v
            # set invers of weight if requested 
            if invw: w = 1. / w

        # take absolut values if requested
            if absolutew: w     = abs(w)

        if absolute:  quant = abs(quant)

        # log of the x axis
        if logx:  
            # check if there are negative values and switch to symlog
                symlog_switch = "symlog"

                # compute ranges and relative absolute values og them 
                aminn         = abs(quant).min()
                amaxx         = abs(quant).max()

                # bins on both sides of zero 
                posbins       =  10**np.linspace(np.log10(aminn), np.log10(quant.max()), bins)
                if (quant<0).any():
                    negbins       = -10**np.linspace( np.log10(-quant.min()),np.log10(aminn), bins)

                    # merge the bin arrays 
                    bins          = np.append(negbins, posbins)
                else:
                    bins          = posbins
                    symlog_switch = "log"
            ##else:
               #

                # regular log bins
                #bins          = 10**np.linspace(np.log10(quant.min()),
                #                      np.log10(quant.max()), bins)
        # plot the histogram 
        w = None
        ax.hist(quant, bins = bins, weights=w, **kwargs)

        # plot additional simple statistics
        if stat:
            mean, std = utils.weighted_avg_and_std(quant, w)
            #ax.axvline(mean,             alpha = salpha    , label="Mean: %.4E"%mean        , color = scolor)
            #ax.axvspan(mean-std,mean+std,alpha = salpha*0.5, label=r"std: $\pm$"+ "%.4E"%std, color = scolor)
            
        if logx:
            # set the scale according to the switch 
            
            ax.set_xscale(symlog_switch)
        # if a label is supplied or statistic is computed, show the legend.
        if "label" in kwargs:
            label = kwargs["label"]
            if label != "" or stat is True:
                ax.legend()
        # set the axes labels 
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        # thigh layout 
        fig.tight_layout()
        return [ax, mean, std] if stat else ax


 ###stuff    
    def q_mean(self, quant:str, typefield:str = "gas", typefieldw:str = "gas", weight:str = "mass" ):
        """
            This function compute the average of the distribution of the quantity q.
            quant       : str, name of the quantity 
            typefield   : str, type of the quantity (gas/DM/nbody...)
            weight      : str, name of the weighting quantity 
            typeweight  : str, name of the weighting quantity 
            return: unyt quantity 
        """
        
        if self.sp==None: self.init_sp()
        self.q= self.sp.quantities.weighted_average_quantity((typefield ,quant),
                                    (typefield, weight))
        return self.q
            
#### DEVIATION 
    def q_std(self, quant:str, typefield:str = "gas", typefieldw:str = "gas", weight:str = "mass" ):
        """
            This function compute the standard deviation of the distribution of the quantity q.
            quant       : str, name of the quantity 
            typefield   : str, type of the quantity (gas/DM/nbody...)
            weight      : str, name of the weighting quantity 
            typeweight  : str, name of the weighting quantity 
            return: unyt quantity 
        """
        if self.sp==None: self.init_sp()
        self.q_dev = self.sp.quantities.weighted_standard_deviation((typefield ,quant),
                              (typefieldw, weight))[0]
        return self.q_dev
    
#### TOTAL QUANTITES
    def eth(self, ):
        if self.sp==None: self.init_sp()
        self.thermal_eng = self.sp.quantities.total_quantity(("gas" ,"thermal_energy")).to("Be")
        return self.thermal_eng
    
    def ek(self, ):
        if self.sp==None: self.init_sp()
        self.kinetic_eng = self.sp.quantities.total_quantity(("gas" ,"kinetic_energy")).to("Be")
        return self.kinetic_eng
    
    def etot(self, ):
        if self.sp==None: self.init_sp()
        th  = self.sp.quantities.total_quantity(("gas" ,"thermal_energy")).to("Be") if self.thermal_eng==None else self.thermal_eng
        kin = self.sp.quantities.total_quantity(("gas" ,"kinetic_energy")).to("Be") if self.kinetic_eng==None else self.kinetic_eng
        self.etot = th+kin
        return self.etot
    
    def radmom(self, ):
        if self.sp==None: self.init_sp()
        self.radial_mom = self.sp.quantities.total_quantity(("gas" ,"radial_momentum")).to("Msun*km/s")
        return self.radial_mom
    def totmass(self, ):
        if self.sp==None: self.init_sp()
        self.radial_mom = self.sp.quantities.total_quantity(("gas" ,"mass")).to("Msun")
        return self.radial_mom
    def full_analysis(self):
        self.thermal_eng=s1=self.eth()    if self.thermal_eng==None else self.thermal_eng 
        self.kinetic_eng=s2=self.ek()     if self.kinetic_eng==None else self.kinetic_eng 
        self.tot_eng=s3=self.etot()       if self.tot_eng==None else self.tot_eng 
        self.radial_mom=s4=self.radmom()  if self.radial_mom==None else self.radial_mom 
        return np.array([s1,s2,s3,s4])
    # YT UTILS
    
    
    def init_sp(self):
        """
            This function initialise the object to directly acces the data with tuples, as in all_data() from yt.
            The reason this is not automatically called is that it can require computational time to load the data, 
            although it is not always needed.
        """
        # if there are no selection parameters, just load all the data
        if self.select_data is None:

            self.sp = self.ds.all_data()
            
        elif self.select_data[0]=="out": 
            # load data outside the given range of values for the given quantity 

            self.data_restr = self.select_data[0]
            self.sp = self.ds.all_data().include_outside(self.select_data[1],
                                                         self.select_data[2],
                                                         self.select_data[3])
        elif self.select_data[0]=="in":
            # same of above, but takes the values inside 

            self.data_restr = self.select_data[0]
            self.sp = self.ds.all_data().include_inside(self.select_data[1],
                                                        self.select_data[2],
                                                        self.select_data[3])  
        else:
            # in any other case
            raise Exception("%r is not allowed."%self.select_data[0] )        
        return

    def add_field_Bratio(self, rho0 = 1., Bc = 2.e-6):
        
        d_units = unyt.gram/unyt.cm**3
        rhoc    = rho0#*d_units
        Bc      = 2e-6 *unyt.Gauss
        def ratio(field, data):
                return (
                    data["gas","magnetic_field_magnitude"].in_units("gauss")/(Bc*(data["gas","density"].d/rho0)**0.5)
                )
        self.ds.add_field(
                name=("gas", "expected_Bratio"),
                function=ratio,
                sampling_type="local",
                display_name=r"$B/[B_0(z_0)(\rho/\rho_0(z_0))^{\alpha}]$",
                units="dimensionless",#"kelvin/cm",
            )
        return
    
    
    def coord(self):
        """
            This function extract the 3D coordinates of each cell. 
            return: list(unyt.array*3), coordinates
        """
        # extract cells coordinates 
        if self.sp == None: self.init_sp()
        if self.x == None:
            self.x = self.sp[("gas","x")]
            self.y = self.sp[("gas","y")]
            self.z = self.sp[("gas","z")]
        return [self.x, self.y,self.z]
    
    def get_cube(self, level:int, field, typefield="gas" , ghost = 0):
        """
            This function extracts a regular grid from the data, at a resolution 2**level.
            level    : int, set the resolution as 2**level
            field    : str, select the field to extract
            typefield: str, select the type of the field
            ghost    : int, how many ghost cells should be used, if needed.
            return: list(regularised 3D quantity, regularized 3D cube for further field extraction) 
        """
        manyfields = False
        if isinstance(field, list): manyfields = True
        if self.arepo: raise Exception("Not available to arepo simulations.")
        # check that the refinement level selected is present in the simulation
        if level < self.base_level: 
            print("l is %d"%level)
            print("The coarser level is l = %d, setting variable level = %d"%(self.base_level,self.base_level))
            level = self.base_level
        elif level > self.base_level + self.maxlevel:
            level = self.base_level+self.maxlevel
            print("l is %d"%level)
            print("The maximum level is l = %d, setting variable level = %d"%((self.base_level + self.maxlevel),
                                                                              (self.base_level + self.maxlevel)))
        
                            
        ref       = 2**level                 # compute the the number of cells per direction
        low       = self.ds.domain_left_edge # set left edge of the simulation box
        dims      = ref                      # redundant but it's ok 
        # use the smoothed_covering_grid function from the yt project. 

        cube      = self.ds.smoothed_covering_grid(level-self.base_level,
                    left_edge = low, dims = dims,
                    fields = [(typefield,field)] if not manyfields else field,
                    num_ghost_zones = ghost)
        u = (cube[field])

        return [u , cube]
    
    
    
    @staticmethod
    def findpos(sp, poss, part):

        """
            This function find the position of the stars in the domain and give then back ordered by index.
            sp    : data container as defined in the class
            poss  : data 
            return: list(np.array*3)
        """
        #if self.sp==None: self.init_sp()
        idx=sp[('nbody',"particle_index")]
        p = np.argsort(idx)
        
        x = poss[0]#*w0
        y = poss[1]#*w0
        z = poss[2]#*w0
        x = x[p]
        y = y[p]
        z = z[p]
        x = x[part]
        y = y[part]
        z = z[part]
        positions = [x,y,z]
        return positions
    
def calc_spectrum_E(k,Kk, kdir, kbins, N):
        whichbin = np.digitize(kdir.flat, kbins)
        ncount = np.bincount(whichbin)
        E = np.zeros(len(ncount) - 1)
        for n in range(1, len(ncount)):
            E[n - 1] = np.sum(Kk.flat[whichbin == n])
        k = 0.5 * (kbins[0 : N - 1] + kbins[1:N])
        k_dir=k[0:len(k)-1]
        E = E[1:N]
        E = E [0: len(E)-1]
        index = np.argmax(E)
        kmax = k[index]
        E_max = E[index]
        return k_dir, E, kmax, E_max
    
def fft_comp(ds, cube,iu, interface, workers, filter, f):
    #rho = cube["gas", "density"].d
    u = cube[iu].d
    nx, ny, nz = u.shape
    if filter!=None:
        #f= cube[filter[1]].d
        if filter[0]=="in":
            
            u[ np.logical_and( f>=filter[2],f<=filter[3])]=filter[4]
            

        elif filter[0]=="out":
            u[ f<=filter[2]]=filter[4]
            u[ f>=filter[3]]=filter[4]
            
        else:
            sys.exit("Only possible keys are 'in' and 'out', entered: %s"%filter[4] )
            
    if interface==scipy:
        
        ru = interface.fft.fftn( u, workers = workers)[
            0 : nx // 2 + 1, 0 : ny // 2 + 1, 0 : nz // 2 + 1
        ]
    else: 
        ru = interface.fft.fftn( u)[
        0 : nx // 2 + 1, 0 : ny // 2 + 1, 0 : nz // 2 + 1
    ]
    norm=1
    norm*=2/nx 
    norm*=2/ny 
    norm*=2/nz 
    ru =   ru *norm
    return np.abs(ru) ** 2



def get_trajectories(   box=1,
                        nstars=1,
                        scale_l   = costs.pc,
                        scale_t   = costs.myr,
                        negative_time = True,
                        row_wise  = False,
                        #dense     = False,
                        data_path = None,
                        base_names= "snp",
                        fmt       = "%02d",
                        ext       = ".dat", 
                        velocities= False,
                        time=True,
                        plot=True):
    
        if isinstance(nstars, int): 
            starlist = range(0, nstars)
        elif isinstance(nstars, list|tuple|np.ndarray ):
            starlist = nstars 
        
        if data_path is None: 
            data_path = os.getcwd()+"/trajectories/"
        else:
            data_path +="/"  
        # initialise grid
        tpath = data_path +  base_names + fmt%(int(nstars)-1) + ext
        with open(file=tpath, mode="r") as f:
                data = np.loadtxt(f)
                ndata = len(data[:,0])
                          #star   time axes   
        stars = np.zeros((nstars,
                          3 if not velocities else 6,
                          ndata)) 
        #print(stars)
        for idx in starlist:
            tpath = data_path +  base_names + fmt%(int(idx)) + ext 
            
            with open(file=tpath, mode="r") as f:
                data = np.loadtxt(f)
            t = data[:,0] #* costs.myr / scale_t
            
            if negative_time: t = t - t[0] 
            x = data[:,1] * costs.pc / scale_l + box * 0.5
            y = data[:,2] * costs.pc / scale_l + box * 0.5
            z = data[:,3] * costs.pc / scale_l + box * 0.5
            if velocities: 
                vx = data[:,4]
                vy = data[:,5]
                vz = data[:,6]
                
            current_t = np.inf
            if row_wise: 
                points_set = [[xi,yi,zi] for (xi,yi,zi) in zip(x[t<= current_t],
                                                   y[t<= current_t],
                                                   z[t<= current_t]) ]
            else:
                
                stars[idx,0,:] = x
                stars[idx,1,:] = y
                stars[idx,2,:] = z
                if velocities:
                    stars[idx,0+3,:] = vx
                    stars[idx,1+3,:] = vy
                    stars[idx,2+3,:] = vz
                   
        
        if not row_wise: points_set=stars
                          
        if plot: 
            fig= plt.figure()
            ax=plt.subplot(projection="3d")
            for star in stars:
                #print(star[0])
                x,y,z = star
                ax.plot(x,y,z)
            #plt.show()
        return points_set  if not time else [points_set, t] 
    
    
    
    
class analyse_RAMSES_clumps(analyse):
    def __init__(self, outpath = "./",
                 base_filename = "output_",
                 fmt           = "%05d", 
                 extention     = "",
                 outnumb:int = -1,
                 #time_series:bool = False, 
                 mhd:bool = False,
                 poisson=False,
                 select_data = None,
                 more_fields = False,
                 #smoothing_factor = 3. # used only for arepo when data are being loaded in the class
                 **kwargs
                 ):
        super().__init__(
                 outpath = "./",
                 base_filename = "output_",
                 fmt           = "%05d", 
                 extention     = "",
                 outnumb = -1,
                 mhd = False,
                 poisson=False,
                 select_data = None,
                 more_fields = False,
                 **kwargs
                 )
        
    def get_ncpu(self):
        infos_file = self.filename + "/info_%05d.txt"%self.outnumb
        with open(infos_file) as f: 
            pass
