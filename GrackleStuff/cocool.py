import numpy as np 
import constants as cc 
import code_units as cd
from yt import load as yt_load
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import root
import unyt
import matplotlib.pyplot as plt
import os
from scipy.optimize import root_scalar
import paper_settings as pg
import logging
import matplotlib.colors as colors
from labellines import labelLine, labelLines
from matplotlib.ticker import ScalarFormatter, NullFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def conductivity_fact(T): return 5.6e-7
paper = pg.PaperSettings()
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=paper.niceColors[1:])



class CoolingFunction:
     
    def __init__(self, path, datapath, write = False):
        self.cached = {}
        self.path = path 
        self.datapath = datapath
        self.false = True
        os.makedirs(path, exist_ok = True)
        self._load_data()
        return
     
    def __load_cached(self):
        compute_cool     = False
        compute_cooltime = False
        compute_rho      = False
        compute_temp     = False
       
        try:
            with open(self.path+"rho_cooling.txt", mode="r") as frho:
                rho = np.array(np.loadtxt(frho, dtype=float))
            self.__cache_update("rho",rho)
            
                #rho = np.array(rho)
        except Exception as e :
            print(e)
            compute_rho = True
            
        try: 
            with open(self.path+"temp_cooling.txt", mode="r") as ftemp:
                temp = np.array(np.loadtxt(ftemp, dtype=float))
            self.__cache_update("temp",temp) 
        except Exception as e: 
            print(e)
            compute_temp = True

        try:
            with open(self.path+"cool.txt", mode="r") as frate:
               rate = np.array(np.loadtxt(frate, dtype=float))
            self.cached.update({"cooling": rate})
        except Exception as e:
            print(e)
            compute_cool = True

        try:
            with open(self.path+"cool_time.txt", mode="r") as fratecool:
                cooltime = np.array(np.loadtxt(fratecool, dtype=float))

            self.__cache_update("cool_time", cooltime)
        except Exception as e:
            print(e)
            compute_cooltime = True

        self.compute = dict(rho  = compute_rho,
                            temp = compute_temp,
                            cool = compute_cool,
                            cool_time=compute_cooltime)
        return
    def __cache_update(self, name, data):
        self.cached.update({name:data})
        return



    def _load_data(self):

        self.__load_cached()
        datasets = os.listdir(self.datapath)
        def tonum(string): return int(string.replace("cooling_rate_","").replace(".h5","") )
        
        datasets.sort(key= lambda name : tonum(name))
        
        compute_rho = self.compute["rho"]
        compute_temp =self.compute["temp"]
         
        compute_cool = self.compute["cool"]
        compute_cooltime = self.compute["cool_time"]
        
        if compute_rho:
            rho = []

            with open(path+"rho_cooling.txt", mode="w") as frho:
                for filename in datasets:
                        print(filename, " COOLING")
                        ds = yt_load(cooling_datasets+filename)
                        d = ds.data["density"][0]/(cc.mu_mol*cc.mH*unyt.g)
                        rho.append(d)
                        frho.write("%.8e\n"%d)                
            rho = np.array(rho)
            self.__cache_update("rho",rho) 
        if compute_cool or compute_temp or compute_cooltime:
            os.system("rm %scool.txt %scool_time.txt %stemp_cooling.txt"%(path,path,path))
            for  filename in datasets:
                    ds = yt_load(cooling_datasets+filename)
                    print(filename,"tmep")
                    #break
                    if first:      
                        temp=ds.data["temperature"]
                        ftemp = open(path+"temp_cooling.txt", mode="+a")
                        for T in temp: 
                            ftemp.write("%.8e\n"%T)
                        tempg,rhog   = np.meshgrid(temp,rho, indexing="xy")
                        self.__cache_update("cooling",np.empty((rho.size,temp.size)))
                        self.__cache_update("cooling_time", np.empty((rho.size,temp.size)))
                        first = False
                    cool      = ds.data["cooling_rate"]
                    # self.__cache_update("cool",cool) 
                    cool_time = ds.data["cooling_time"]
                    cool_time = cool_time/abs(cool_time)
                    # self.__cache_update("cool_time",cool_time) 
                    #fcool.write(cool)
                    with open(path+"cool.txt", mode="ab") as frate:
                        np.savetxt(frate, cool.reshape(1, cool.shape[0]), fmt='%.8e')
                    with open(path+"cool_time.txt", mode="ab") as fratetime:
                        np.savetxt(fratetime, cool_time.reshape(1, cool_time.shape[0]), fmt='%.8e')
                    self.cooling[i]     = cool
                    self.cooling_time[i]= cool_time
        return
    @property
    def cooling(self): return self.cached["cooling"]
    @property
    def rho(self): return self.cached["rho"]
    @property
    def temp(self): return self.cached["temp"]
    @property
    def cooltime(self): return self.cached["cool_time"]

path = "/work/pacicco/analysis_playground/GrackleStuff/mycooling_datasets/HM/"
cooling_datasets = path+"/cooling_datasets/"

ccool = CoolingFunction(path=path, datapath = cooling_datasets)
cooling = ccool.cooling 
cooling_time = ccool.cooltime 
rho = ccool.rho
temp = ccool.temp

#print(rho.min(), rho.max(), temp.min()

tmax = temp.max()
first = True
equcoord = [] 
#fcool = open("rho_cooling.txt", mode="+a")

def find_eq(rho, temp, cooltime, i):
    qplus  = np.where(cooltime>0)[0]#np.argmin(cooltime)
    qminus = np.where(cooltime<0)[0]#np.argmin(cooltime)
    print(qplus, qminus)
    #print(rho[i], temp[qminus],temp[qplus], cooling[qminus], cooling[qplus])
    #print(rho[i], temp[qmin], cool)            
            
            


logrho,logtemp   = np.log10(rho), np.log10(temp)
logcooling       = np.log10(cooling)
tempg,rhog       = np.meshgrid(temp,rho, indexing="xy")
logtempg,logrhog = np.meshgrid(logtemp,logrho, indexing="xy")


interpLogLambda_nT   = RegularGridInterpolator(points= (logrho,logtemp),values = logcooling)
interpLogLambda_time = RegularGridInterpolator(points= (logrho,logtemp),values = cooling_time)


def calcT(n,P):
    return P/(cc.kB*n)
def calcP(n,T):
    return n*cc.kB*T
def calcn(P,T): return P/(cc.kB*T)    

def LogLambda_nT(logn, logT):
    return  interpLogLambda_nT((logn, logT))

def LogLambdaTime_nT(logn, logT):
    return  interpLogLambda_time((logn, logT))


def Lambda_nT(n, T):
    #print(n.min(),n.max(), T.min(), T.max()/1e8)
    T[T>tmax] = tmax
    logn, logT = np.log10(n), np.log10(T)
    return 10**LogLambda_nT(logn, logT)

def LambdaTime_nT(n, T):
    logn, logT = np.log10(n), np.log10(T)
    return LogLambdaTime_nT(logn, logT)




if __name__=="__main__":
    

        
    RT_eq = []

    for ni in rho:
    
        def wrap(TT):
            return LambdaTime_nT(ni,TT)
        try:
            tt = root_scalar(wrap, bracket=[temp.min(), temp.max()], x0=1e4,
                            #xtol=1e-8,
                            #rtol=1e-8, 
                            maxiter=1000)
            RT_eq+=[[ni,tt.root]]
        except:
            RT_eq+=[[-1, -1]]
            
        pass
    RT_eq   =  np.array(RT_eq)
    q   = np.where(RT_eq>-1)[0]
    RT_eq = RT_eq[q]
    n_eq    =  RT_eq[:,0] #10**RT_eq[:,0] 
    T_eq    =  RT_eq[:,1] #10**RT_eq[:,1]
    P_eq     =  calcP(n_eq, T_eq)

    Logn_eq = np.log10(n_eq)
    LogT_eq = np.log10(T_eq)
    LogP_eq = np.log10(P_eq)

    G = cc.factG_in_cgs

    def rhof(n):
        return cc.mu_mol * cc.mH * n  

    def cooling_time(n, T,):
        gamma = 5./3. 
        P = calcP(n, T)
        return  cc.kB*T/( (gamma-1)*n * Lambda_nT(n,T) )   

    #    return  (1./(gamma-1))*P/( n * Lambda_nT(n,T) )   

    def speed_of_sound(T ):
        gamma = 5./3. 
        return np.sqrt(gamma *cc.kB * T / (cc.mu_mol * cc.mH)) 

    def cooling_length(n,T ):

        return cooling_time(n, T )*speed_of_sound(T)

    def jeans_length(n, T):
        num = np.pi * speed_of_sound(T)**2
        den = G * rhof(n)
        return np.sqrt(num / den)

    def free_fall_time( n ):
        return np.sqrt((3.0*np.pi)/(32.0 * G * rhof(n)))

    def conductivity(T):
        return conductivity_fact(T) * T**(5./2.)


    def conduction_time(n,T, L, gamma = 5./3.):
        P = calcP(n, T)
        return 7./(2*(gamma-1.0))*P/(conductivity(T)*T/L**2)


    def field_length(n, T):
        num = conductivity(T) * T 
        den = (n**2 * Lambda_nT(n,T))
        return np.sqrt(num/den)




    plt.rcParams.update({
    "text.usetex": True, **paper.mplfigfontstyle
    #"font.family": "Helvetica"
    })
    #plt.rcParams.update({'font.size': 9})

    n0 = 7.005273279999999492e-01 
    T0 = 6.167396659348973117e+03
    P0 = calcP (n0, T0) 
    nh = rho
    T  = temp
    netCool = 10**LogLambda_nT(logrhog,logtempg).T 
    figxxx = plt.figure(12)
    figxxx.set_size_inches( paper.onecol,paper.onecol)
    axxx = figxxx.add_subplot(111)
    h    = axxx.pcolormesh( nh,T, netCool,cmap="BuGn"
                           ,linewidth=0,rasterized=True,
                           norm = colors.LogNorm())
 
    axxx.plot(n_eq,T_eq, color = "black", lw = paper.figure.linewidth,
               linestyle = "dashed",
               label = r"$n_\mathrm{eq},T_\mathrm{eq}$",
               alpha = 0.7
              )
    axxx.legend(loc = "lower left")

    cbaxes = inset_axes(axxx, width="70%", height="5%", loc="upper center") 
    
    cb = figxxx.colorbar(h, cax = cbaxes,
                    orientation="horizontal",
                    label=r"$\log_{10}{(\|\Lambda(n,T)\|)}$",
                    extend='both')
    # set colorbar label plus label color
    cb.set_label(label=r"${\|\Lambda(n,T)\|}$",color="white")
    # set colorbar tick color
    cb.ax.yaxis.set_tick_params(color="white")
    cb.ax.xaxis.set_tick_params(color="white")
    # set colorbar edgecolor 
    cb.outline.set_edgecolor("white")

    # set colorbar ticklabels
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color="white")
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color="white")
    axxx.set_xlabel(r"$n_{\rm H}\,\,[{\rm cm^{-3}}]$")
    axxx.set_ylabel(r"$T\,\,[{\rm K}]$")
    
    
    axxx.set_xscale("log")
    axxx.set_yscale("log")
    figxxx.tight_layout()
    figxxx.savefig(path+"Grackle_full_cool.pdf" ,dpi  = paper.figure.dpi )# dpi = np.rint(np.ceil(nfiles/paper.onecol)))
    figxxx.savefig(path+"Grackle_full_cool.png" , dpi = paper.figure.dpi)
    figxxx.savefig(path+"Grackle_full_cool.jpeg", dpi = paper.figure.dpi)
    fig = plt.figure(0)
    fig.set_size_inches( paper.onecol+0.2*paper.onecol,paper.twocol)
    
    ax0 = fig.add_subplot(111)
    ax0.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False)
    ax0.set_xticks([])
    
    h   = ax0.pcolormesh(nh, T, netCool,cmap="BuGn", norm=colors.LogNorm())

    ax0.plot(n_eq, T_eq, color = "black", lw=2, label = r"$\rho_{\rm eq},T_{\rm eq}$ curve")
    
    #fig.colorbar(c, ax = ax0, label=r"$\log_{10}(\|\Lambda-\Gamma\|) [\rm erg\,\rm cm^{3}$]")
    #cc.kB=1
    ax0.legend(loc = "lower left")
    #fig.tight_layout()
    cbaxes = inset_axes(ax0, width="70%", height="5%", loc="upper center") 
    
    fig.colorbar(h, cax = cbaxes, orientation="horizontal", label=r"$\log_{10}{(\|\Lambda(n,T)\|)}$")#,location="top")
        
      
    #ax0.scatter(n0, T0,       marker = "x", color = "red", zorder = 1000, label="Initial conditions")
    ax0.set_ylabel(r"$T\,\,[\rm K]$")
    
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    
    
    figg = plt.figure(1)
    figg.set_size_inches( paper.onecol,paper.onecol)
    ax1 = figg.add_subplot(111)
    ax1.plot(n_eq, P_eq/cc.kB, color = "black", lw=2, label = r"$\rho_{\rm eq},T_{\rm eq}$")
    ax1.fill_between(n_eq, P_eq/cc.kB, color = "red", alpha=0.1) 
    ax1.fill_between(n_eq, P_eq/cc.kB,y2=1.05*(P_eq/cc.kB).max(), color = "blue", alpha=0.1)       
    with open(path+"equilibrium_NT.dat", mode="w") as f:
        #for nn,tt in zip(n_eq, T_eq): f.write("%.15f %.15f \n"%(nn, tt) )
        np.savetxt(f, np.column_stack((n_eq, T_eq)))
    
    with open(path+"equilibrium_NP.dat", mode="w") as f:
        #for nn,pp in zip(n_eq, P_eq): f.write("%.15f %.15f \n"%(nn, pp) )
        np.savetxt(f, np.column_stack((n_eq, P_eq)))
    
    with open(path+"equilibrium_RHOp.dat", mode="w") as f:
        #for nn,tt in zip(n_eq*cc.mu_mol*cc.mH, P_eq): #f.write("%.15f %.15f \n"%(nn, tt) )
        np.savetxt(f, np.column_stack((n_eq*cc.mu_mol*cc.mH, P_eq)))
    
    #ax1.scatter(n0, P0/cc.kB, marker = "x", color = "red", zorder = 1000, label="Initial conditions")
    ax1.set_ylabel(r"$P\,\,[\rm P/k_b]$")
    
    ax1.set_xlabel(r"$n_{\rm H}\,\,[\rm cm^{-3}]$")
    pmin, pmax = (P_eq/cc.kB).min(), (P_eq/cc.kB).max()#ax1.get_ylim()
    nmin, nmax = (n_eq).min(), (n_eq).max()#ax1.get_ylim()
    
    
    labels_xvalues = []
    rhott = rho  
    for texp in range (1, 5 ):
        tt = 10**texp
        
        Pp = calcP(rhott, tt) /cc.kB   
        pp0 = 1500.0#10**((np.log10(pmin)+np.log10(pmax))*0.2) 
        #print(calcn(pp0*cc.kB, tt), tt)
        labels_xvalues+=[calcn(pp0*cc.kB, tt)]
        tmprho = rhott #np.log10(rhott)      #[Pp>=pmin]
        #tmprho = tmprho[Pp>=pmin]
        #Pp     = Pp    [Pp>=pmin]
        #tmprho = tmprho[Pp<=pmax]
        #Pp     = Pp    [Pp<=pmax]

        ax1.plot(tmprho, Pp, label = r"$%i$ K"%tt, zorder=1, lw = 2)
    ax1.set_ylim(0.95*pmin, 1.05*pmax)
    ax1.set_xlim(0.95*nmin, 1.0*nmax)
    

    labelLines(ax1.get_lines()[1:], xvals=labels_xvalues, zorder=2.5)
    #labelLines(ax1.get_lines(), zorder=2.5)

    
    #ax1.legend(ncols = 2)

    

    ax1.set_xscale("log")
    ax1.set_yscale("log")
     
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    ax1.xaxis.set_minor_formatter(ScalarFormatter())
    
    figg.tight_layout()
    figg.savefig(path+"Pressure-Density.pdf", dpi=512)
    figg.savefig(path+"Pressure-Density.png", dpi=512)
    figg.savefig(path+"Pressure-Density.jpeg",dpi=512)
    fig.subplots_adjust(hspace=0.)
    fig.savefig(path+"grackcool.png", dpi = paper.figure.dpi)
    quit()
    fig2 = plt.figure(2)
    ax21   = fig2.add_subplot(221)
    ax22   = fig2.add_subplot(223)    
    ax211  = fig2.add_subplot(222)
    ax222  = fig2.add_subplot(224)
    
    coolrate_eq = Lambda_nT(n_eq, T_eq) 
    tc    = cooling_time   (n_eq, T_eq)
    lc    = cooling_length (n_eq, T_eq)
    lJ    = jeans_length   (n_eq, T_eq)
    tff   = free_fall_time (n_eq)
    lF    = field_length   (n_eq, T_eq)
    tcond = conduction_time(n_eq, T_eq, lF)
    
    tc = tc/cc.kyr
    tff = tff/cc.kyr
    tcond = tcond / cc.kyr
    
    lc = lc/cc.pc
    lJ = lJ/cc.pc
    lF = lF/cc.pc
    
    ax21.plot(n_eq, tc   , label = r"$\tau_c$"         , color = "crimson"     )
    ax21.plot(n_eq, tff  , label = r"$\tau_{\rm tff}$"  , color = "seagreen"    )
    ax21.plot(n_eq, tcond, label = r"$\tau_{\rm cond}$", color = "midnightblue")
    
    ax22.plot(n_eq, lc   , label = r"$\lambda_c$"      , color = "crimson"    )
    ax22.plot(n_eq, lJ   , label = r"$\lambda_{\rm J}$", color = "seagreen"    )
    ax22.plot(n_eq, lF   , label = r"$\lambda_{\rm F}$", color = "midnightblue")
    
    ax211.plot(T_eq, tc   , color = "crimson"     )
    ax211.plot(T_eq, tff  , color = "seagreen"    )
    ax211.plot(T_eq, tcond, color = "midnightblue")
    
    ax222.plot(T_eq, lc   , color = "crimson"     )
    ax222.plot(T_eq, lJ   , color = "seagreen"    )
    ax222.plot(T_eq, lF   , color = "midnightblue")
    
    
    
    ax21.set_ylabel(r"$\tau\,\,[\rm {kyr}]$")
    ax22.set_ylabel(r"$\lambda\,\,[\rm {pc}]$")
    ax22.set_xlabel(r"$n_{\rm H}\,\,[{\rm cm^{-3}}]$")
    ax222.set_xlabel(r"$T\,\,[{\rm K}]$")
    ax21.legend()
    
    ax21.axvline( n0, label = r"$n_0 = %.1f\,\,{\rm cm^{-3}}$"%n0, color="black", linestyle="dashed", lw = 3, alpha = 0.5)
    ax22.axvline( n0, label = r"$n_0 = %.1f\,\,{\rm cm^{-3}}$"%n0, color="black", linestyle="dashed", lw = 3, alpha = 0.5)
    ax211.axvline(T0, label = r"$T_0 = %.1f\,\,{\rm K}$"%T0, color="black", linestyle="dashed", lw = 3, alpha = 0.5)
    ax222.axvline(T0, label = r"$T_0 = %.1f\,\,{\rm K}$"%T0, color="black", linestyle="dashed", lw = 3, alpha = 0.5)
    

    fig2.subplots_adjust(hspace = 0.1 , wspace = 0.1 )   
    
    ax22.legend()
    ax222.legend()
    ax21.set_xscale("log")
    ax22.set_xscale("log")
    ax21.set_yscale("log")
    ax22.set_yscale("log")
    ax211.set_xscale("log")
    ax222.set_xscale("log")
    ax211.set_yscale("log")
    ax222.set_yscale("log")
    
    
    ax21.set_xticks([])
    ax211.set_yticks([])
    ax222.set_yticks([])
    ax211.set_xticks([])
    def inistudd(n0,T0):       
        coolrate0 = Lambda_nT(n0, T0)
        tc0    = cooling_time   (n0, T0)
        lc0    = cooling_length (n0, T0)
        lJ0    = jeans_length   (n0, T0)
        tff0   = free_fall_time (n0)
        lF0    = field_length   (n0, T0)
        tcond0 = conduction_time(n0, T0, lF0)
    
        tc0    /= (1e5*cc.kyr)
        lc0    /= (1e5*cc.kyr)
        tcond0 /= (1e5*cc.kyr)
        lJ0    /= cc.pc
        tff0   /= cc.pc
        lF0    /= cc.pc
        
        
        
        print("Cooling rate %.2e"%coolrate0)
        print("N0, T0 = %.2f, %.2f"%(n0, T0))
        print("tc    %.4f 10^5 kyr  "%(tc0))
        print("tff   %.4f 10^5 kyr  "%(tff0))
        print("tcond %.4e 10^5 kyr  "%(tcond0))

        print("lc    %.4f pc"%(lc0))
        print("lJ    %.4f pc"%(lJ0))
        print("lF    %.4f pc"%(lF0))
        
    
    inistudd(n0,T0)
#    print("")
#    f = 0.1
#    n0t = n0*(1.-f)
#    T0 = calcT(n0t, P0)
#    inistudd(n0t,T0)
#    print("")
#    n0t = n0*(1.+f)
#    T0 = calcT(n0t, P0)
#    inistudd(n0t,T0)
#    print("")
#    f = 0.9
#    n0t = n0*(1.-f)
#    T0 = calcT(n0t, P0)
#    inistudd(n0t,T0)
#    print("")
#    n0t = n0*(1.+f)
#    T0 = calcT(n0t, P0)
#    inistudd(n0t,T0)
    
    plt.show()

    












