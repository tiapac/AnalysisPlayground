import numpy as np

def close_factors(number):
    ''' 
    find the closest pair of factors for a given number
    '''
    factor1 = 0
    factor2 = number
    while factor1 +1 <= factor2:
        factor1 += 1
        if number % factor1 == 0:
            factor2 = number // factor1
        
    return factor1, factor2
def almost_factors(number):
    '''
    find a pair of factors that are close enough for a number that is close enough
    '''
    while True:
        factor1, factor2 = close_factors(number)
        if 1/2 * factor1 <= factor2: # the fraction in this line can be adjusted to change the threshold aspect ratio
            break
        number += 1
    return factor1, factor2
def compute_rows_columns(npart, tipo, orientation):
    if tipo=="square":
        numb= almost_factors(npart)
    elif tipo=="rectangle":
        numb = almost_factors(npart)
        if orientation=="vertical":
            numb= max(numb), min(numb)
        elif orientation=="horizontal": 
            numb= min(numb), max(numb)
    return numb 
def friediagram( cs, ca, plot=True,ax=None, fig=None, lw=5, numpint=500):
        import matplotlib as plt

        cms2=ca**2 + cs**2
        theta=np.linspace(0,np.pi, numpint )
        Vfast=np.sqrt(0.5 * (cms2 + (cms2**2 - 4*ca**2*cs**2*(np.cos(theta))**2)**0.5))
        Vslow=np.sqrt(0.5 * (cms2 - (cms2**2 - 4*ca**2*cs**2*(np.cos(theta))**2)**0.5))
        Vint=ca*abs(np.cos(theta))
        d=[theta, Vfast, Vint, Vslow]
        if plot: 
            if fig == None: fig=plt.figure()
            if ax  == None: ax=fig.add_subplot(projection="polar")
            ax.set_thetamin(0) # set the limits
            ax.set_thetamax(180) 
            ax.set_xlabel(r"$v\,\,\rm{[km/s]}$")
            ax.set_title(r"$\theta\,\,\left[\frac{v \cdot B}{ \|v\|\|B\|} \right]$", pad=-50)
            #ax.xaxis.set_label_coords(0.5, 0.075) # change the location of the xlabel to given x, y locations w.r.t. the entire figure
            ax.plot(theta, Vfast, color="midnightblue", label="Fast")
            ax.plot(theta, Vint, color="green", label="Intermidiate")
            ax.plot(theta, Vslow, color="red",  label="Slow")
            ax.plot(theta, [ca]*numpint, linestyle="--", color="black", label=r"$C_A$")
            ax.plot(theta, [cs]*numpint, linestyle= ":", color="black", label=r"$C_S$")
            ax.plot(theta, [cms2**0.5]*numpint,linestyle="-.", color="black", label=r"$(C_A+C_S)^{1/2}$")
            ax.legend()
            result= [ax,d,fig] 
            
        else:
            result=d
        return result