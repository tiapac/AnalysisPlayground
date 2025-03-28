import numpy as np
import os
from utils.utils import *
import time as Time
from paper_settings import PaperSettings
import time as TIME


paper = PaperSettings()
paper.timeStampLoc = (0.1,0.3)
myunits = paper.units
ave   = True
path  = "./"

savepath = path+"/my_elliaxy2/"
if not os.path.exists(savepath): os.mkdir(savepath)

start    = 1
end      = 400
plot     = True
analyse  = True
namebase = "axes_length_"



projection  = [4, 3]
xyz         = "xyz"
axes_string = "yzyx"
axes        = { pj:axes_string[pj-1] for pj in projection}

xyz_ind     = {1:"x",2:"y",3:"z"}
xyz_to_ind  = {"x":1,"y":2,"z":3}

axes_pairs  = [find_axes(axxx) for axxx in  axes_string ]

kind        = ["dens","dens", "dens"]
cmaps       = [paper.cmaps.colormaps[k] for k in kind]

log         = [True]*len(kind)
nbins       = [50]*len(kind) 
scales      = {"dens":"scale_d",
               "temp":"unit",
               "|v|" :"scale_b",
               "press":"unit",
               "pmag":"unit"}
clim  = [[1e-27,5e-22],[1e-27,5e-22],[1e-27,5e-22],[1e2,1e8]]#[[1e-28, 1e-22],[None, None],[1e-28, 1e-22],[None, None]]#*len(kind) 

#for key in aliases.keys():
#    myunits.labels.update({key:myunits.labels[aliases[key]]})



update_always = True
initime   = Time.time()


def checktime(msg=""): print("--- %s %.2e seconds ---" %(msg,(Time.time() - initime)))  

checktime("Start")

if plot:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib import collections  as mc
    plt.rcParams.update({
                "text.usetex": True, **paper.mplfigfontstyle
                #"font.family": "Helvetica"
                })
    nplots = len(projection)
    if nplots%2!=0:
        nrows = 1
        ncols = nplots
    else: 
        nrows = nplots//2
        ncols = 2 if nplots >1 else 1
    fig, axs = plt.subplots(nrows = nrows, ncols = ncols , **{"figsize":(paper.onecol * ncols, paper.onecol * nrows)})
    axs = axs.flatten()# axs.flatten()
    fig.set_size_inches(paper.twocol, paper.twocol)
else: 
    axs = range(len(kind))
if analyse: 
    filey = open(file=savepath+namebase+"y.dat", mode="w")
    filex = open(file=savepath+namebase+"x.dat", mode="w")
    filez = open(file=savepath+namebase+"z.dat", mode="w")
    for num in range(start, end, 1):

        for ind, ax in enumerate(axs):
            current_axis = axes[projection[ind]]
            try:
                if plot: ax.clear()
                data, [time, Lx, Ly, Lz],[frame_nx, frame_ny], info = load_map(path,  proj = projection[ind], kind = kind[ind], num = num )
                
                density, _, _, _ = [data, None, None, None] if kind[ind] == "dens" else load_map(path,  proj = projection[ind], kind = "dens", num = num )
                temp   , _, _, _ = [data, None, None, None] if kind[ind] == "temp" else load_map(path,  proj = projection[ind], kind = "temp", num = num )
                        
                dx = max( [Lx, Ly, Lz] ) / frame_nx
                dy = max( [Lx, Ly, Lz] ) / frame_ny
                cd       = code_scales(info)
                data     *= cd[scales[kind[ind]]]
                old_info = dict(info)
                boxlen   = info["boxlen"]
                
                xedges, yedges   = np.arange(1, frame_nx+1), np.arange(1, frame_ny+1)
                xc  = (xedges - 0.5) * dx - 0.5 * boxlen
                yc  = (yedges - 0.5) * dx - 0.5 * boxlen
                
                x, y    = np.meshgrid(xc, yc, indexing = "xy")
                        
                m = density*dx**3
                        
                data_flat = data.flatten()    
                x_flat    = x.flatten() 
                y_flat    = y.flatten()  
                m_flat    = m.flatten()
                temp_flat = temp.flatten()
                
                xmin, xmax, ymin, ymax = x_flat.min(), x_flat.max(), y_flat.min(), y_flat.max()
                
                cond      = temp_flat > 1e5
                data_flat = data_flat[cond]
                x_flat   = x_flat[cond]
                y_flat   = y_flat[cond]
                z_flat = np.zeros_like(x_flat) 

                m_flat = m_flat[cond]
                temp_flat = temp_flat[cond]
                
                M  = np.sum(m_flat)
                xm = np.sum((m_flat*x_flat))/M
                ym = np.sum((m_flat*y_flat))/M
                zm = np.sum((m_flat*z_flat))/M

                dxm = x_flat - xm
                dym = y_flat - ym
                dzm = z_flat - zm
                
                I_xx = np.sum(m_flat * (dym**2 + dzm**2))
                I_yy = np.sum(m_flat * (dxm**2 + dzm**2))
                I_zz = np.sum(m_flat * (dxm**2 + dym**2))
            
                I_xy = -np.sum(m_flat * dxm * dym)
                I_xz = -np.sum(m_flat * dxm * dzm)
                I_yz = -np.sum(m_flat * dym * dzm)

                I = np.array([[I_xx, I_xy, I_xz],
                              [I_xy, I_yy, I_yz],
                              [I_xz, I_yz, I_zz]])    
                
                eigenvalues, eigenvectors = np.linalg.eigh(I)
                
                a_squared = (5.0 / (2.0 * M) ) * (eigenvalues[1] + eigenvalues[2] - eigenvalues[0])
                b_squared = (5.0 / (2.0 * M) ) * (eigenvalues[0] + eigenvalues[2] - eigenvalues[1])
                c_squared = (5.0 / (2.0 * M) ) * (eigenvalues[0] + eigenvalues[1] - eigenvalues[2])
                axes_lengths = np.sqrt([a_squared, b_squared, c_squared]) #np.sqrt([a_squared, b_squared, c_squared])
                a, b, c = axes_lengths#*0.9
                if plot:
                    cff = ax.imshow(data, cmap = cmaps[ind],
                                extent=(xmin,xmax,ymin,ymax ),origin="lower", norm = LogNorm())#
                
                    ax.scatter(xm, ym, s = 10, color= "red")
                p = np.argsort(eigenvalues)
                eigenvalues = eigenvalues[p]
                eigenvectors = eigenvectors[:, p]*1  # Sort eigenvectors accordingly
                semi_axes = [
                    (xm + a * eigenvectors[0, 0], ym + a * eigenvectors[1, 0]),  # Major axis (aligned with eigenvector 0)
                    (xm + b * eigenvectors[0, 1], ym + b * eigenvectors[1, 1]),  # Minor axis (aligned with eigenvector 1)
                ]

                # Prepare the lines to plot
                lines = [
                    [(xm, ym), semi_axes[0]],  # Major axis
                    [(xm, ym), semi_axes[1]],  # Minor axis
                ]
                if plot:
                    lc = mc.LineCollection(lines, colors = "red", linewidths = 2)
                    ax.add_collection(lc)
                print(axes[projection[ind]], ind)
                if current_axis == "y":
                    ordered = (num, time, (a+b)/2. )
                    filey.write("%d %f %f \n"%ordered)#(projection[ind], num, time, a,b,c))
                
                elif current_axis == "x":
                    
                    ordered = (num, time, min(a,b), max(a,b))
                    filex.write("%d %f %f %f \n"%ordered)#(projection[ind], num, time, a,b,c))
                elif current_axis == "z":
                    ordered = (num, time, min(a,b), max(a,b))
                    filez.write("%d %f %f %f \n"%ordered)#(projection[ind], num, time, a,b,c))
                
                #print(f"radius done for {num:05d}")
            
                if plot: fig.savefig(savepath+"axis_%05d.png"%num)
            except: 
                continue
    filex.close()
    filey.close()
    filez.close()

    
for ax in "xyz": 
    with open(file = savepath+namebase+"%s.dat"%ax, mode ="r") as f: 
        types = "int,float,float" if ax=="y" else "int,float,float,float"
        data = np.loadtxt(f,dtype = types)
        
           
        
       
#fig.savefig("averad.png")



