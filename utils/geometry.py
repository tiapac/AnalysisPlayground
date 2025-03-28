import numpy as np 



class coordinates:
    def __init__(self, *args):
        self.ndims = ndims = len(args)
        self.x = args[0]
        if ndims>1:
            self.y = args[1]
            if ndims>2: self.z = args[2]
        pass
    def filter_coordinates(self,*args,**kwargs):
        if self.ndims == 2 :
            return self.filter_coordinates2D(*args,**kwargs)
        else:
            return NotImplementedError(f"Filter_coordinates works only in Only 2D. ndims =  {self.ndims}")

    def filter_coordinates2D(self, theta, inverted = False):
        """
        Filters points within the angular range -theta to theta with respect to the x-axis.

        Parameters:
            x (numpy.ndarray): x-coordinates of the grid.
            y (numpy.ndarray): y-coordinates of the grid.
            theta (float): Maximum angle (in radians) from the x-axis.

        Returns:
            numpy.ndarray: Masked array of indices of the filtered points.
        """
        x, y = (self.x, self.y) if not inverted else (self.y, self.x)

        # Compute angles
        phi = np.arctan2(y, x)

        # Create a mask for angles within the range
        mask1 = (phi >= -theta) & (phi <= theta)

        # Mask for specular range
        mask2 = (phi >= np.pi - theta) & (phi <= np.pi + theta)
        mask3 = (phi >= -np.pi - theta) & (phi <= -np.pi + theta)

        # Combine masks
        mask = mask1 | mask2 | mask3

        # Return filtered coordinates
        return mask


