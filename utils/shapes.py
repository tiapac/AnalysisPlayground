import numpy as np 

class shape:

    def __init__(self):
        pass
    def drawEllipsoid(X, SemiAxes, ngrid = 50):
        xCenter, yCenter, zCenter = X
        a, b, c = SemiAxes
        #draw sphere
        u, v = np.mgrid[0:2*np.pi:ngrid, 0:np.pi:ngrid]
        x = np.cos(u)*np.sin(v)
        y = np.sin(u)*np.sin(v)
        z = np.cos(v)
        # shift and scale sphere
        x = a * x + xCenter
        y = b * y + yCenter
        z = c * z + zCenter
        return (x,y,z)
    def drawSphere(X, r, ngrid = 50):
        return shape.drawEllipsoid(X, (r,r,r), ngrid = ngrid)