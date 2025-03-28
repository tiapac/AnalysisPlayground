import matplotlib.pyplot as plt
import os 
import sys
from scipy.interpolate import CubicSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import interpolate
from scipy import stats
import numpy as np 
import pandas as pd
import warnings

class stellar_tracks:

    def __init__(self, 
              metals_path = "./evo_traks/tables_grids2011/Z014/",
              rotation    = True,
              start       = None,
              end         = None,
              outfilename = "./STELLAR_EVO_TABLE_polyfit_new_test.txt"
            ):
        """
            metals_path : str, path to the data.
            rotation    : bool, choose or not the data for the rotating start
            start/end   : int|float, choose the initial/final mass range for the fit.
                                    ints are taken as indexes while floats as masses 
                                    used to pick the index of the simulation with closest
                                    mass 
            outfilename: str       , name of the filenames where to print the results.
        """





        # initialise names of the files and related masses
        self.metals_path = metals_path
        self.filenames, self.masses = self.set_names(path = metals_path, rotation = rotation)
        if start is None: start = 0
        if end is None:   end   = len(self.filenames)-1
        # set initial and final indexes 
        if isinstance(start, float): self.start = np.argmax(self.masses >= start)#list(self.masses).index(start)
        if isinstance(end, float): 
            self.end = np.argmax(self.masses >= end)
            if self.end > len(self.filenames): self.end = len(self.filenames) 

        self.outfilename = outfilename
        if os.path.exists(self.outfilename): os.system("rm "+self.outfilename)
        
        self.names     = ["mass", "lg(Teff)","lg(L)","Gamma_Ed","lg(Md)"]
        self.ram_names = ["mass","teff","lum","gamma_ed","mdot"]
        self.labels =[r"$M_\bigstar$", r"$\rm Log(T_{\rm eff})$", r"$\rm Log(L_\bigstar/L_\odot)$", r"$\Gamma_{\rm ed}$", r"$M_{\rm loss}$"]


    def main_loop(self, plot = True):
        """
            Loop over the various files and extract the data.
            SEEMS OK BUT NEEDS TESTING
        """
        figs = [] if plot else None
        for indx0,filename in enumerate(self.filenames):
            if indx0 >= self.start and indx0 <= self.end:
                if plot:
                    fig = plt.figure(indx0, figsize=(8,8))
                    ax  = fig.add_subplot(111)
                # import data
                data = pd.read_csv(self.metals_path+filename, header=0, skiprows=[1,2],sep="\s+")# delim_whitespace=True)
                time = np.float64(data["time"])
                with open(self.outfilename, "a") as fi:
                    self.ini_sub_text(indx0, fi, self.masses, self.start, self.end)
                    masssplines=[]
                    teffsplines=[]
                    lumsplines= []
                    gammasplines=[]
                    mdotsplines=[]
                    arrays_list = [masssplines,teffsplines,lumsplines,gammasplines,mdotsplines]
                    for (indx,i),name, array in zip(enumerate(self.names),self.ram_names, arrays_list):    
                        data_tofit = data[i]
                        u, c = np.unique(time, return_index=True)    
                        
                        warnings.simplefilter('ignore', np.RankWarning)
                        #print(u,c)
                        

                        for j in range(time.size-1):
                            
                            knots=np.array([time[j],time[j+1]])
                            knotsy=np.array([data_tofit[j],data_tofit[j+1]])
                            
                            if i == "mass": knotsy = np.log10(knotsy)
                            m, q = np.polyfit(knots, knotsy   ,deg=1)
                            
                            base = ""
                            ends  = ""      
                        
                            if indx==0:
                                if j==0 :
                                    base  = "  IF (tt .LE. %+.9f) THEN \n "%(knots[1])
                                elif j < time.size-2:
                                    base  = "  ELSE IF ( tt .GT.%+.9f .AND. tt .LE. %+.9f) THEN \n "%(knots[0],knots[1])

                                else:
                                    base  = "  ELSE IF ( tt .GT.%+.9f) THEN \n "%(knots[1])

                            if indx==4 and j==time.size-2:
                                    if indx0 <= self.end: 
                                        ends = "\n  ENDIF"
                                    if indx0 == self.end: ends +=" \nENDIF\n\nendsubroutine choose_model"
                                            

                                    
                            array += self.data_to_append(j, name, array, time, knotsy, q, m, indx, base, ends )
                                
                        if plot:
                            s=4
                            if i!="mass" and i!= "Gamma_Ed":
                                ax.plot(time/1e6, data_tofit, label=self.labels[indx], lw=s)
                                gg=1
                            elif  i== "Gamma_Ed":
                                    ax.plot(time/1e6,  np.log10(data_tofit), label=self.labels[indx], lw=s)
                            elif i=="mass":
                                #gg=1
                            #if i=="mass":
                                ax.plot(time/1e6,  np.log10(data_tofit), label=self.labels[indx], lw=s)
                                
                                fig.suptitle(r"Model for mass $%f$" %self.masses[indx0])

                            ax.legend()
                            ax.set_xlabel(r"Time [Myr]")
                            #fig.savefig("Plot_model_%f.png"%self.masses[indx0])
                    data_to_write = np.column_stack([masssplines,teffsplines, lumsplines, gammasplines, mdotsplines])
                    np.savetxt(fi , data_to_write, fmt="%s")
                    # print closure
                figs += [fig]
        return figs
                    
    @staticmethod
    def data_to_append(j, name, array, time, knotsy, q, m, indx, base, ends):
        """
            This function append data to the arrays.
        """
        comma = ";" if (indx != 1 and indx != 4)  else "\n"
        if j==0 :
            string = ( "%s=%+.9f %s" %(name,knotsy[0], comma))
        elif j == time.size-2:
            string = ( "%s=%+.9f %s" %(name,knotsy[1], comma))
        else:
            string = ( "%s=%+.9f %+.9f*tt %s" %(name, q, m, comma))

#        array+=
        return [base + "  "+ string + ends]

    @staticmethod
    def set_names(path, rotation):
        """
            Set names of the files needed for the fit and the corresponding masses.
        """
        if rotation is True: rotation = 4 
        lista = os.listdir(path)
        mass_ini = []
        namelist = []
        for indx,name in enumerate(lista):
            lista[indx] = name.replace("p", ".")
            name_ = lista[indx]
            rotation = int(name_[len(name_)-5])
            if (rotation == 4 ):##only with rotation

                namelist += [name]
                mass_ = float(name_[1:4])
                mass_ini += [mass_]
        # sort 
        mass_ini = np.array(mass_ini)
        namelist = np.array(namelist)
        p        = np.argsort(mass_ini)
        mass_ini = mass_ini[p]
        namelist = namelist[p]

        return namelist, mass_ini
    @staticmethod 
    def ini_sub_text(i0, file, masses, start, end):
        """
            Print main stuff in the subroutine body.
        """
        if i0 == start:
                        file.write("subroutine choose_model(m_model, tt, mass, teff, lum,mdot, gamma_ed)\n\
    use amr_commons, only: dp\n\
    implicit none\n\
    real(dp):: tt\n\
    integer::m_model\n\
    real(dp)::mass, teff, lum, mdot, gamma_ed\n\n\
IF ( m_model == %i  ) THEN !!for mass %f\n"%(i0+1, masses[i0]) )
        elif i0 > start and i0 <= end:
            file.write("ELSE IF ( m_model == %i  ) THEN !!for mass %f\n"%(i0+1, masses[i0]) )
        return

if __name__ == "__main__": 
    st = stellar_tracks(start = 0., end = 120.)
    st.main_loop(plot=True)
    plt.show()
                    