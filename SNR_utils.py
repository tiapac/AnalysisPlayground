import constants as cc 
import numpy as np 
from utils.shapes import shape

drawEllispoid = shape.drawEllipsoid
drawSphere = shape.drawSphere
ESN = 1e51
E51 = ESN/1e51

XH  = 0.76
msg = "Not allowed to use this time, sholud be smaller then the starting of the phase."



def rho_n(n):
    return n*cc.mu_mol*cc.mH
    
def Rst(t, n, dummyxi = None):
    """
    Sedov phase raius, in parsec
    t_: time in yr
    rho_0: denisty of the gas in g cm^-3 

    """
    rho=rho_n(n)
    
    xi5=2.026
    return (xi5*1.e51/rho)**(1./5.)*(cc.yr*t)**(2./5.)/cc.pc




def tsed(n, xi =None, mej=15, vej=None):
    mej *= cc.M_sun
    vej = np.sqrt(2.0 * ESN / mej) if vej is None else vej * cc.kms
    rho = rho_n(n)
    R3 = (3.0 * mej) / (rho * 4.0 * np.pi)
    R = R3**(1.0 / 3.0)
    t = R / vej
    return t / cc.yr
    
def t_pds(n, xi):
    """
    n: density of the gas in cm^-3
    xi: metallicity of the gas in solar
    """
    nH      = XH*n
    t_sf=3.61e4*(E51**(3./14.))/(nH**(4./7.)*(0.132*xi/0.218)**(5./14.))
    return t_sf/np.e

def Rpds(t , n, xi):

    tpds = t_pds(n,xi)
    if t<tpds:  print(UserWarning(msg))
    return Rst(tpds , n)*(t/tpds)**(2./7.)

def t_mcd(n,xi,v_1e4_kms = None ):
    nH      = XH*n
    if v_1e4_kms is None: v_1e4_kms = 1.0 
    phi_spitzer= 1.0
    t1 = 61.0*v_1e4_kms**3/((0.132*xi/0.218)**(9./14.)*nH**(3./7.)*E51**(3./14.))
    t2 = 476./(0.132*xi/0.218*phi_spitzer)**(9./14.)

    return min(t1,t2)*t_pds(nH,xi)

def Rmcd(t,n,xi):
    nH      = XH*n
    tmcd = t_mcd(n,xi)
    if t<tmcd:  print(UserWarning(msg)) 
   
    return Rpds(tmcd,n,xi)*(t/tmcd)**(1./4.)

def t_merge(n,xi, cs5 = 10.0):
    nH      = XH*n
    beta = 2.0
    cs6 = cs5/10.0
    return 153.*((E51**(1./14)*nH**(1./7.)*\
    (0.132*xi/0.218)**(3./14.)/(beta*cs6)))**(10./7.)*\
    t_pds(n, xi)

def Rmerge(t, n, xi , cs5 = 10.0):
    nH      = XH*n
    tmerge  =  t_merge(n,xi, cs5)
    if t<tmerge:  print(UserWarning(msg)) 
    tpds    =  t_pds(n,xi)
    tmcd    =  t_mcd(n,xi)
    if tmerge <= tpds:
        radius =   Rst
    elif tmerge > tpds and tmerge <= tmcd :
        radius =   Rpds
    elif  tmerge > tmcd:
        radius =   Rmcd 
    return radius(tmerge,n,xi)+max(cs5*cc.kms*(t-tmerge)*cc.yr/cc.pc,0.0)

def composit_R(t, n, xi , cs5 = 10.0):
    
    tpds    =  t_pds(n,xi)
    tmcd    =  t_mcd(n,xi)
    tmerge  =  t_merge(n,xi, cs5)
    #print(tpds, tmcd, tmerge, t)
    if t>tmerge:
        radius = Rmerge(t, n, xi,cs5)
    else:
        if t <= tpds:
            radius =   Rst(t,  n, xi)
        elif t > tpds and t <= tmcd :
            radius =   Rpds(t, n, xi)
        elif  t > tmcd:
            radius =   Rmcd(t, n, xi)
    return radius

def RadiusSNR(t, n, xi , cs5 = 10.0):

    return np.array([composit_R(ti, n, xi, cs5 = cs5) for ti in t])

def composit_colors(t, n, xi , cs5 = 10.0):
    tpds    =  t_pds(n,xi)
    tmcd    =  t_mcd(n,xi)
    tmerge  =  t_merge(n,xi, cs5)
    if t>tmerge:
        radius = 4
    else:
        if t <= tpds:
            radius =   1
        elif t > tpds and t <= tmcd :
            radius =   2
        elif  t > tmcd:
            radius =   3
    return radius        


def RSB_Badry(t,n,theta = 0.3,deltat_SNe = 5e5):
    # t is in yr, converted to Myr, deltaSNe in SN is in Myr and converted to 0.1Myr
        
    return 83*((1.0-theta)*E51)**(1./5.)*((deltat_SNe/0.1e6)*n)**(-1./5.)*(t/1e6)**(3./5.)
def RSB_Fede(t,n,NSN, deltaOB=4.e7):
    t=t*cc.yr
    nH  = n*XH
    #Lsn = 3.17e36*NSN*E51/(deltaOB/1.0e7)
    Lsn = NSN*1e51/(20*cc.myr)
    
    # t is in yr, converted to Myr, deltaSNe in SN is in Myr and converted to 0.1Myr
    
    #rsb = 1.67e2*(Lsn/1e37)**(1./5.)#*(t/1e7)**(3./5.)
    #rsb = rsb /(nH**(-1./5.))
    #rsb = rsb * (t/1e7)**(3./5.)
    rsb = 0.8*(Lsn/rho_n(n))**(1./5.)*t**(3./5.)/cc.pc
    return rsb 
if __name__=="__main__":
    x=RSB_Badry(1.2e6,0.429,theta = 0.0,deltat_SNe = 1e6)
    x1=RSB_Fede(1.2e6,0.429,2)
    #rsb=
    print(x1,x)
    