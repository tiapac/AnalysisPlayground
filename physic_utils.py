import numpy as np 
import constants as cc 
import code_units as cd
from yt import load as yt_load
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import root
import unyt
import matplotlib.pyplot as plt
import os
import logging
#ds = yt_load("./cooling_rateUVB.h5")       
#temp = ds.data["temperature"]
#rho  = ds.data["density"]/(cc.mu_mol*cc.mH*unyt.g)
inifiles=0
nfiles = 800

write= False
compute_cool = False
compute_cooltime = False
compute_rho  = False
compute_temp = False

path = "/work/pacicco/analysis_playground/cooling_extracted/"
cooling_datasets = "/work/pacicco/analysis_playground/cooling_datasets_zoom/"
try:
    os.makedirs(path)
except: 
    pass

try:
    with open(path+"rho_cooling.txt", mode="r") as frho:
        rho = np.loadtxt(frho, dtype=float) 
    rho = np.array(rho)
except Exception as e :
    print(e)
    compute_rho = True
try: 
    with open(path+"temp_cooling.txt", mode="r") as ftemp:
                    temp = np.loadtxt(ftemp, dtype=float) 
    temp=np.array(temp) 
except Exception as e: 
    print(e)
    compute_temp = True
try:
    with open(path+"cool.txt", mode="r") as frate:
                    cooling = np.loadtxt(frate, dtype=float)
except Exception as e:
    print(e)
    compute_cool = True
try:
    with open(path+"cool_time.txt", mode="r") as fratecool:
                    cooling_time = np.loadtxt(fratecool, dtype=float)
except Exception as e:
    print(e)
    compute_cooltime = True

if compute_rho or write:
    rho = []

    with open(path+"rho_cooling.txt", mode="w") as frho:
        for i in range(inifiles,nfiles):
                print(i)
                ds = yt_load(cooling_datasets+"cooling_rate_%02d.h5"%i)
                d = ds.data["density"][0]/(cc.mu_mol*cc.mH*unyt.g)
                rho.append(d)
                frho.write("%.8e\n"%d)                
    rho = np.array(rho)
first = True
#fcool = open("rho_cooling.txt", mode="+a")
if compute_cool or compute_temp or compute_cooltime or write:
    os.system("rm %scool.txt %scool_time.txt %stemp_cooling.txt"%(path,path,path))
    for i in range(inifiles,nfiles):
            ds = yt_load(cooling_datasets+"cooling_rate_%02d.h5"%i)
            print(i)
            #break
            if first:      
                temp=ds.data["temperature"]
                ftemp = open(path+"temp_cooling.txt", mode="+a")
                for T in temp: 
                    ftemp.write("%.8e\n"%T)
                tempg,rhog=np.meshgrid(temp,rho, indexing="xy")
                cooling=np.empty((rho.size,temp.size))
                cooling_time=np.empty((rho.size,temp.size))
                first = False
            cool = ds.data["cooling_rate"]
            cool_time = ds.data["cooling_time"]
            cool_time = cool_time/abs(cool_time)
            #fcool.write(cool)
            with open(path+"cool.txt", mode="ab") as frate:
                np.savetxt(frate, cool.reshape(1, cool.shape[0]), fmt='%.8e')
            with open(path+"cool_time.txt", mode="ab") as fratetime:
                np.savetxt(fratetime, cool_time.reshape(1, cool_time.shape[0]), fmt='%.8e')
            cooling[i]=cool
            cooling_time[i]=cool_time
            
            
def calcT(n,P):
    return P/(cc.kB*n)
def calcP(n,T):
    return n*cc.kB*T
          

logrho,logtemp = np.log10(rho), np.log10(temp)
logcooling     =  np.log10(cooling)
tempg,rhog=np.meshgrid(temp,rho, indexing="xy")
logtempg,logrhog=np.meshgrid(logtemp,logrho, indexing="xy")
Pg = calcP(rhog,tempg)
logPg= np.log10(Pg)

interpLogLambda_nT   = RegularGridInterpolator(points= (logrho,logtemp),values = logcooling)
interpLogLambda_time = RegularGridInterpolator(points= (logrho,logtemp),values = cooling_time)



def LogLambda_nT(logn, logT):
    return  interpLogLambda_nT((logn, logT))

def LogLambdaTime_nT(logn, logT):
    return  interpLogLambda_time((logn, logT))


def Lambda_nT(n, T):
    logn, logT = np.log10(n), np.log10(T)
    return 10**LogLambda_nT(logn, logT)

def LambdaTime_nT(n, T):
    logn, logT = np.log10(n), np.log10(T)
    return LogLambdaTime_nT(logn, logT)

    

    
G = cc.factG_in_cgs

def rho(n):
    return cc.mu_mol * cc.mH * n  

def cooling_time(n, T, cooling_rate ):
    
    P_n = cc.kB*T

    tc = P_n/( n * cooling_rate )   
    return  tc

def speed_of_sound(T, gamma = 5./3. ):
    return np.sqrt(gamma *cc.kB * T / (cc.mu_mol * cc.mH)) 

def cooling_length(n,T, cooling_rate):

    return cooling_time(n, T, cooling_rate )*speed_of_sound(T)

def jeans_length(n,T):
    num = np.pi * speed_of_sound(T)**2
    den = G * rho(n)
    return np.sqrt(num / den)

def free_fall_time( n ):
    return np.sqrt((3.0*np.pi)/(32.0 * G * rho(n)))

def conductivity(T):
    return 5.6e-7 * T**(5./2.)

def field_length(n, T, cooling_rate):
    num = conductivity(T) * T 
    den = (n**2 * cooling_rate)
    return np.sqrt(num/den)

if __name__=="__main__":
    n0 = 1e-2
    T0 = 3e3
    coolrate0 = Lambda_nT(n0,T0) 
    tc = cooling_time(n0, T0, coolrate0)
    lc = cooling_length(n0, T0, coolrate0)
    js = jeans_length(n0,T0)
    ff = free_fall_time( n0 )
    cs = speed_of_sound( T0)
    lf = field_length(n0, T0,coolrate0)
    tc = tc/cc.kyr
    lc = lc/cc.pc
    js = js/cc.pc
    ff = ff/cc.myr
    cs = cs/cc.kms
    lf = lf/cc.pc
    
    print("cool rate = %.2e myr"%coolrate0 )
    print("tc = %.2e kyr"%tc )
    print("lc = %.2e pc "%lc )
    print("js = %.2e pc "%js )
    print("ff = %.2e myr"%ff )

    print("cs = %.2e kms"%cs )
    print("lf = %.2e pc "%lf )
    fig = plt.figure(0) 
    ax  = fig.add_subplot(111)
    if 1: coolint=ax.pcolormesh(           
                logrhog,
                logtempg,
                LambdaTime_nT(rhog,tempg),
                #LogLambda_nT((logrhog,logtempg)),
                #vmin=-30,
                #cmap="kamae"
                )
    plt.show()
    quit()
    
    plot_cT= LambdaTime_nT(rhog,calcT(rhog, Pg))

    
    coolint=ax.pcolormesh(           
                logrhog,
                np.log10(10**logPg/cc.kB),
                plot_cT,
    #            LogLambda_nT((logrhog,logtempg)),
                #vmin=-30,
                cmap="kamae"
            )
    rho0 = 1.0*(cc.mu_mol*cc.mH)
    T0   = 0.3e4
    P0   = rho0/(cc.mu_mol*cc.mH)*cc.kB*T0 
    
    def Tt(var)  :
       r1 = rho0*(1.0-var)
       r2 = rho0*(1.0+var)
       
       return  P0/(r1/(cc.mu_mol*cc.mH)*cc.kB), P0/(r2/(cc.mu_mol*cc.mH)*cc.kB) 
     
    def RHOt(var): 

        return rho0*(1-var), rho0*(1+var)
    
    #ax.scatter(np.log10(rho0)*(cc.mu_mol*cc.mH)        ,np.log10(P0/cc.kB))
    #ax.scatter(np.log10(RHOt(0.1)[0]/(cc.mu_mol*cc.mH)),np.log10(P0/cc.kB))
    #ax.scatter(np.log10(RHOt(0.1)[1]/(cc.mu_mol*cc.mH)),np.log10(P0/cc.kB))

    #fig = plt.figure(0) 
    #ax  = fig.subplot(111)
    
    
    pass