try:
    import analysis_playground.constants as cc
    #from analysis_playground.constants import * 
except:
    class __cc:
        gamma        = 5./3.
        factG_in_cgs = 6.6740800e-8
        pc    = 3.0856776e18
        kpc   = 1e3 * pc 
        kpc2cm = kpc
        yr    = 3.15576000e+07
        myr   = yr * 1e6
        M_sun = 1.9891e33
        kms   = 1e5
        mu_mol       = 1.2195e0
        mH           = 1.6605390e-24
        kB           = 1.3806490e-16
        twopi        = 6.2831853e0
        pi           = twopi/2e0

        kB_SI        = 1.3806490e-23
        e_SI         = 1.60217663e-19
    cc = __cc()

from utils.logger import setup_logger

logger = setup_logger("clump_structures")
try:
    import matplotlib.pyplot as plt
    #plt.style.use('dark_background')

    plt_ok = True

except Exception as e:
    logger.warning("%r ||| No matplotlib package found."%e)
    plt_ok = False
try:
    import pyvista as pv
except Exception as e:
    logger.warning("%r ||| No pyvista package found."%e)
       
try:
    # apparently math.sqrt is faster.
    import math
    sqrt = math.sqrt 
except:
    def sqrt(value): 
        return value**(0.5) 
try: 
    import operator
except ImportError: 
    sort_key = lambda x: x.ID # use a lambda if no operator module
else: 
    sort_key = operator.attrgetter("ID") # use operator since it's faster than lambda
    
import numpy as np
from   pandas import read_csv





logger.info("Starting clump analyser.")

def check_axes(ax, fig):
    if not plt_ok:
            logger.error("Cannot plot without matplotlib. It needs to be installed.")
            raise ModuleNotFoundError
    if fig is None: 
            existing = plt.get_fignums()
            if len(existing)>1: 
                num = max(existing)+1
            elif len(existing)==1:
                num = existing[0]+1
            else:
                num = 0
            fig = plt.figure(num)
    if ax is None: 
            ax  = fig.add_subplot(111, projection="3d")
    return ax, fig

 

     
class leaf: 
    """
        Store information related to a leaf. 
    """
    def __init__(self,
                ID,                
                x, y, z,
                rho   = None,
                T     = None, 
                emag  = None,
                ekin  = None,
                ink   = None,
                level = None,
                simulation= None
                ) -> None:
        self.ID        = ID
        self.x           = x 
        self.y           = y 
        self.z           = z
        self.level       = level 
        self.rho    =   rho
        self.T      =   T*cc.mu_mol
        self.emag   =   emag
        self.ekin   =   ekin
        self.ink    =   ink
        self.simulation       = simulation
        self.m = self.mass()
        self.Bmag = self.B()
        return
 
    def mass(self):
        boxsize = (self.simulation.boxlen*self.simulation.code_units.scale_l)
        dx = boxsize/(2**self.level)
        vol= dx**3 
        m = self.rho*vol 
        return m
    def B(self):
        
        emag = self.emag/(self.simulation.code_units.scale_v**2.) #8*np.sqrt(self.emag/self.simulation.code_units.scale_v**2.)*self.simulation.code_units.scale_b
        b = sqrt(2*emag)*self.simulation.code_units.scale_b
        return b
 
 
 
 
    
class clump:
    """
        Store informations related to a clump. 
    """
    def __init__(self,
                ID, 
                x         = None,
                y         = None,
                z         = None,
                level     = None,
                parent    = None,
                ncells    = None,                
                rho_minus = None,
                rho_plus  = None,
                rho_av    = None,
                mass_cl   = None,
                relevance = None,
                simulation= None,
                leafs = []) -> None:
        
        self.ID          = ID 
        self.level       = level
        self.parent      = parent
        self.ncells      = ncells
        self.x           = x 
        self.y           = y 
        self.z           = z
        self.rho_minus   = rho_minus 
        self.rho_plus    = rho_plus 
        self.rho_av      = rho_av
        self.mass_cl     = mass_cl 
        self.relevance   = relevance
        self.leafs       = leafs
        self.simulation  = simulation
        #self.nleafs      = len(self.leafs)
        
        # to be defined
        self.r           = None
        pass
    
    def add_leaf(self, __leaf:leaf|list) -> None:
        """
            Ad a leaf or a list of leafs to the clump
        """
        if self.leafs == []:
            if isinstance(__leaf,list|tuple):
                self.leafs = __leaf
            else:
                self.leafs = [__leaf]    
        else:
            if isinstance(__leaf,list|tuple):
                self.leafs.extend(__leaf)
            else:
                self.leafs.append(__leaf)
        
        return 
    def set_simulation(self, simu):
        self.simulation = simu
        return
    def leafs_data(self):
        """
        Retreive coordinates of the leafs
        """
        return [lef.x for lef in self.leafs],[lef.y for lef in self.leafs],[lef.z for lef in self.leafs],[lef.level for lef in self.leafs]
    def map(self, name):
           return [getattr(lef, name) for lef in self.leafs]
    
    def avg(self, name):
        nleafs = 0
        average_value = 0.0
        for lef in self.leafs:
            average_value+= getattr(lef, name)*lef.m      
            nleafs+=1     
        average_value/=self.mass()#nleafs
        return average_value
    
    def jeans_length(self):
        rho = self.avg("rho")
        cs2  = cc.gamma*cc.kB*self.avg("T")/(cc.mH * cc.mu_mol)
        
        return np.sqrt(np.pi * cs2/(rho*cc.factG_in_cgs))
    def mass(self):
        """
            Boxlen is in p
        """
        m = 0.0
        boxsize = (self.simulation.boxlen*self.simulation.code_units.scale_l)
        for lef in self.leafs:
            dx = boxsize/(2**lef.level)
            vol= dx**3 
            m+= lef.rho*vol
        return m 
    def radius(self, X0 = [0.0]*3):
        """
            Compute the distance of the clump from the point X0.
        """
        #if self.r is None: 
        result = sqrt((self.x - X0[0])**2.0 + (self.y - X0[1])**2.0 + (self.z - X0[2])**2.0)
        #else:       
        
        return result
    
    def plot(self,ax = None, fig = None,bounds=None, level = None, map = None,offscreen=True,cbar=False,cm=False, cmap = "viridis",log=True,pyvista=False, **kwargs):
        if not pyvista: 
            ax, fig = check_axes(ax, fig)
        else:
            if fig is None: fig = pv.Plotter(off_screen = offscreen)
            ax = None
        logger.trace("Plotting %i of clump #%i on level %r"%(self.ID, len(self.leafs), "all" if level is None else level))
        x, y, z, lvl = self.leafs_data()
        if cm:
            x0,y0,z0=self.center_of_mass()
            x-=x0
            y-=y0
            z-=z0
            
            
            
        if map is not None: 
            m = np.array(self.map(map))#np.array([clp.ID]*len(x))#np.array([clp.ID]*len(x))
        if level is not None:
            x = np.array(x)
            y = np.array(y)
            z = np.array(z)
            lvl = np.array(lvl)
            x   = x[lvl==level] 
            y   = y[lvl==level] 
            z   = z[lvl==level] 
            m   = m[lvl==level] 
            
            lvl = lvl[lvl==level]
        if log: m = np.log10(m)
        if pyvista:
            points = np.column_stack((x, y, z))
            point_cloud = pv.PolyData(points)
            point_cloud[map] = m
            fig.add_points(
                    point_cloud,
                    scalars   = map,
                    cmap      = cmap,
                    style     = "points_gaussian",#"points", #
                    emissive  = True ,
                    log_scale = log,
                    scalar_bar_args={"color":"white"},
                    **kwargs
                    )
            fig.set_background("black")
            if bounds is None:
                X=1.5*max(abs(x.min()),abs(x.max()))
                Y=1.5*max(abs(y.min()),abs(y.max()))
                Z=1.5*max(abs(z.min()),abs(z.max()))
                L=max(X,Y,Z)
            else:
                L=bounds
            bounds=[-L,L,-L,L,-L,L]
            fig.show_bounds(        bounds= bounds ,
                                    color     ="white",
                                    bold      =False,
                                    location  ='outer',                       
                                    ticks     ='both',                       
                                    n_xlabels =5,                        
                                    n_ylabels =5,                        
                                    n_zlabels =5,                        
                                    xtitle    ="x [pc]",                       
                                    ytitle    ="y [pc]",                      
                                    ztitle    ="z [pc]",    
                                    font_size = 20,
                                    
                                    ##font_size=20                 
                                    )
           
        else:
            if map is None:
                hbar = ax.scatter(x,y,z, **kwargs)#, c = m, cmap="viridis")
            else: 
                hbar = ax.scatter(x,y,z, c = m, cmap=cmap, **kwargs)
                label= "Log10(%s)"%map if log else map
                if cbar: fig.colorbar(hbar, label=label)#,cax=ax)
            if not ax.get_legend_handles_labels() == ([], [], []):
                ax.set_zlabel("z")
                ax.set_xlabel("x")
                ax.set_ylabel("y")
             
        
        return fig, ax
    def center_of_mass(self):#, X0=[0.0,0.0,0.0]):
        #in cgs
        cm = [0.0,0.0,0.0]
        boxsize = (self.simulation.boxlen)*self.simulation.code_units.scale_l
        for lef in self.leafs:
            dx = boxsize/(2**lef.level)
            vol= dx**3 
            m  = lef.rho*vol
            cm[0]+= lef.x*m
            cm[1]+= lef.y*m
            cm[2]+= lef.z*m
        cm[0]=(cm[0]/self.mass())#/self.simulation.units.scale_l
        cm[1]=(cm[1]/self.mass())#/self.simulation.units.scale_l
        cm[2]=(cm[2]/self.mass())#/self.simulation.units.scale_l
        self.x=cm[0]
        self.y=cm[1]
        self.z=cm[2]
        return cm 
    def __len__(self):
        return 0 if self.leafs is None else len(self.leafs)
    
    def __iter__(self):
        """
            Loop over the the lead of a clump.
        """
        return iter(self.leafs)
    
    def __contains__(self, item:leaf):
            
        return item in self.leafs 
    
    def __getitem__(self, index):
         return self.leafs[index]
    
    
    
# just to list them all 
class clumps:
    """
        Just store the clumps. 
    """
    min_cells = 30
    def __init__(self, simulation = None) -> None:
        self.clumps = []
        self.simulation = simulation
        self.leaves     = None
        pass
    def avg(self, name, scale=1.0):
        return [clp.avg(name)/scale for clp in self.clumps]    
    def mass(self, scale=cc.M_sun):
        return [clp.mass()/cc.M_sun for clp in self.clumps]
    def jeans_length(self, scale=cc.pc):
        return [clp.jeans_length()/cc.pc for clp in self.clumps]
           
    
    def read_leafs_data(self):
        logger.info("Reading Leaves data from file...")
        ncells = 0
        ncpus  = self.simulation.ncpu
        outnumb=self.simulation.outnumb
        boxlen =self.simulation.boxlen
        path   = self.simulation.path 
        myleafs  = leafs()
        for cpu in range(1, ncpus+1):
            
            data_path = path+"output_%05d/clump_map.csv%05d"%(outnumb,cpu)
            data = read_csv(data_path,engine="c",dtype=float, names=["x","y","z","level", "ID","rho", "T", "emag","ekin","ink"], header = None)
            id, x, y, z, level, rho, T, emag, ekin, ink = data["ID"],data["x"],data["y"], data["z"], data["level"],data["rho"],data["T"],data["emag"],data["ekin"],data["ink"]
            for i, xx, yy, zz,lvl,rhoo, Tt, emagg, ekinn, inkk in zip(id, x, y, z, level,rho, T, emag, ekin, ink):
                ncells+=1  
                
                myleafs.add(leaf(i,
                                    xx - boxlen/2.0,
                                    yy - boxlen/2.0,
                                    zz - boxlen/2.0,
                                    level = lvl,
                                    rho   = rhoo, 
                                    T     = Tt,
                                    emag  = emagg,
                                    ekin  = ekinn, 
                                    ink   = inkk, 
                                    simulation=self.simulation
                                    ))    
            logger.trace("Reading cpu %i, loaded %i leaves."%(cpu, ncells))
        logger.info("Found %i leaves in total."%len(myleafs))
        
        self.leaves = myleafs
        self.leaves.sort()
        logger.debug("Leaves sorted by index.")
        return self.leaves
            
    def assign_leaves(self):
        logger.info("Creating clumps from leaves.")
        if self.leaves is None: logger.warning("There are no leaves loaded in the clump group.")
        last_leaf_id  = 0
        nclumps       = 0
        # EXPLOIT SORTING
        for n in range(0, len(self.leaves)):
            lef = self.leaves[n] 
            leaf_ID = lef.ID
            if n == 0 or (leaf_ID != last_leaf_id):
                nclumps += 1
                
                self.add(clump(ID = leaf_ID)) 
                self.clumps[nclumps-1].add_leaf(lef)
                logger.trace("Created clump #%i"%leaf_ID) 

            else:    
                self.clumps[nclumps-1].add_leaf(lef)
            last_leaf_id = int(leaf_ID)

        # remove clumps with less then n cells
        ini_length = len(self.clumps) 
        #myclumps[:] = [clp for clp in myclumps if ((len(clp) > min_cells) & (clp.avg("T")>0) )]
        self.clumps[:] = [clp for clp in self.clumps if (len(clp) > self.min_cells)]
        for clp in self.clumps:
            clp.set_simulation(self.simulation)
            clp.center_of_mass()
        logger.info("%i/%i clumpls with less than %i leaves have been removed."%(len(self.clumps), ini_length, self.min_cells))
        return
    def plot(self,fig = None, ax = None, condfun = lambda x: True, **kwags):
        ax, fig = check_axes(ax, fig)
        logger.info("Plotting the clumps in 3D.")
        for clp in self:
            if condfun(clp): clp.plot(fig = fig, ax = ax, **kwags)
            
        return fig, ax
    def plotAclump(self, ID, ax = None, fig = None, **kwargs):
        
        if isinstance(ID, list|tuple):
            logger.debug("Plotting clumps: %r"%ID)
            for clp in self:
                if clp.ID in ID: 
                    logger.trace("Plotting clump #%i"%clp.ID)
                    clp.plot(ax = ax, fig=fig, **kwargs)
        
        else:
            logger.debug("Plotting clump #%i"%ID)
            ax, fig = check_axes(ax, fig)
            for clp in self:
                if clp.ID == ID: clp.plot(ax = ax, fig=fig,  **kwargs)
        return fig, ax
    def set_min_cells(self, value):
        """
            Set the minimum cells without overwriting the property for all instances. 
        """
        self.min_cells = int(value) 
        return    
     
    def add(self, __clump : clump) -> None:
        """
            Add a clump to the list. 
        """
        self.clumps.append(__clump)
        
        return
    def remove(self, __clump):
        
        self.clumps.remove(__clump)
        return
    def sort(self):
        """
            Sort by ID
        """
        self.clumps.sort(key = sort_key)
        return
    
    def __iter__(self):
        """
            Loop over the clump list.
        """
        return iter(self.clumps)
    
    def __len__(self):
        return len(self.clumps)
    
    def __contains__(self, item): 
        """
            Check if the ID of the object is already in the clump list. Works
            for leafs and clumps. Can be slow. 
        """     
        for item in self.clumps:
            if item.ID == self.clumps:
                return True 
        return False 
    def __getitem__(self, index):
        """
            Define indexing. Not related to climps' ID.
        """
        return self.clumps[index]
     
    def __setitem__(self, index, __clump):
        """
            Define index assignement. 
        """
        self.clumps[index] = __clump

class leafs:
    """
        Just to store several leafs.
    """
    
    def __init__(self) -> None:
        
        self.leafs = []
        pass
    
    def add(self,__leaf:leaf):
        self.leafs.append(__leaf)
        return
    
    def sort(self):
        """
            Sort by ID.
        """
        self.leafs.sort(key = sort_key)
        return
            
    def __iter__(self):
        return iter(self.leafs)
    
    
    def __len__(self):
        return len(self.leafs)
    
    
    def __getitem__(self, index):
        """
            To index the leafs
        """
        return self.leafs[index]
   