import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

import MAS_library as MASL
import Pk_library as PKL

from scipy.fft import fftfreq
from scipy.fft import fft
from scipy.fft import ifft
from sys import argv


#################################################################################################################
for i in range(N_iter):
    #np.savetxt("iteration_" + str(i), pos)
    pos = zeldovich_reversed(pos, delta, deltagrid, BoxSize)
    delta = delta_from_pos(pos, deltagrid, BoxSize, n_grid_b_0)
#################################################################################################################
np.savetxt("delta_end_511", delta[511])

posR = np.sqrt((pos[:,0]-BoxSize/2)**2 + (pos[:,1]-BoxSize/2)**2 + (pos[:,2]-BoxSize/2)**2)
posMU = (pos[:,2]-(BoxSize/2))/posR
posPHI = np.arctan2(pos[:,1]-BoxSize/2, pos[:,0]-BoxSize/2)

oldmask2 =  (posR >= minR) & (posR <= maxR) & (posMU >= minMU)

posZ = r_to_z(posR)

np.savetxt(cat_name + "_glass_" + argv[1] + ".dat", np.c_[posZ[oldmask2], posMU[oldmask2], posPHI[oldmask2]])

#oldmask2 =  (posR >= minR) & (posR <= maxR)
#np.savetxt(cat_name + "_glass_full_" + argv[1] + ".dat", np.c_[posZ[oldmask2], posMU[oldmask2], posPHI[oldmask2]])
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################

