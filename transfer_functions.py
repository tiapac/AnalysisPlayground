import numpy as np
from pyvista import opacity_transfer_function as transfer 
def gau(sigma, t, mu):
    
    return np.exp( -0.5 * (t - mu)**2 / sigma**2 )

def tff( model = None, i = 32, ncolors=64, alpha=1.5):
    
    #print(model)
    if model == "sin":    
    
        n    = 3
        t    = np.linspace(0, np.pi * n, ncolors)+np.pi
        func = np.sin(3*t)*(1 - t / t.max())**( alpha)
        func[func<0.] = 0.
        
    elif model == "gaussian":
        
        t     = np.linspace(0, 40, ncolors)
        sigma = 4
        mu    = i #t.max()/i# WORKING FOR LINEAR SCALE              
        func  = gau(sigma, t, mu) / gau(sigma, mu, mu)
    
    elif model is None:
    
        func  = np.ones(1)
    elif isinstance(model, float) or isinstance(model, int):
        func  = np.ones(1)*model
    else:
        # this is the case it is a specific array
        func = model 
        
    
    
    opacity   = transfer(func, 256).astype(float)

    return  opacity
