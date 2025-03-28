import numpy as np
from constants import *
import unyt
scale_l      = 3.0856776e18# * unyt.cm
scale_m      = 1.9885e33   # * unyt.g
scale_d      = scale_m / scale_l**3 #* unyt.g / unyt.cm**3
#scale_t converts time from user units into seconds
scale_t      = 1. / np.sqrt( factG_in_cgs * scale_d ) #* unyt.s
scale_v      = scale_l / scale_t #* unyt.cm/unyt.s
scale_nH =  0.70651e0 / mH * scale_d


scale_b = (scale_t /(np.sqrt(4.* np.pi* scale_d) *scale_l))**(-1)