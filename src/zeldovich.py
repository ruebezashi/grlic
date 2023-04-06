import gc
import numpy as np
from numba import jit

import warnings
from numpy.fft import rfftn, irfftn, fftfreq, rfftfreq



import MAS_library as MASL

from settings import cell_size

#calculates the fourier space displacement field from the fourier space density contrast
@jit(nopython=True)
def psi_from_delta(delta, kx, ky, kz, ncells):
    len_x = delta.shape[0]
    len_y = delta.shape[1]
    len_z = delta.shape[2]
    psi = np.empty((len_x,len_y,len_z,3), dtype = np.complex64)
    for i in range(len_x):
        for j in range(len_y):
            for k in range(len_z):
                if ((i == 0) and (j == 0) and (k == 0)) or ((i == int(ncells[0]/2)) or (j == int(ncells[1]/2)) or k == (int(ncells[2]/2))):
                    psi[i,j,k,:] = 0
                else:
                    psi[i,j,k,:] = - 1j * np.array([kx[i],ky[j],kz[k]]) / (kx[i]**2 + ky[j]**2 + kz[k]**2) * delta[i,j,k]
    return(psi)


def displace_reverse(pos, delta, box_size):
    print("displacing particles...")

    ncells = (box_size/cell_size).astype("int") + (box_size/cell_size).astype("int")%2

    # fourier transform delta in 3D                                                                                                                                                                                                                                                                                           
    delta = rfftn(delta).astype("csingle")

    #unit of fftfreq is in cycles / Mpc/h -> to get to h/Mpc multiply by 2pi                                                                                                                                                                                                                                                  
    kx = (fftfreq(ncells[0], box_size[0]/ncells[0])*2*np.pi).astype("float32")
    ky = (fftfreq(ncells[1], box_size[1]/ncells[1])*2*np.pi).astype("float32")
    kz = (rfftfreq(ncells[2], box_size[2]/ncells[2])*2*np.pi).astype("float32")
    #get FT-displacement field from FT-delta                                                                                                                                                                                                                                                                                  
    
    #reverse displacement field
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        psi = psi_from_delta(delta, kx, ky, kz, ncells)
    del delta
    gc.collect()

    #inverse FT of displacement Field                                                                                                                                                                                                                                                                                         
    psi = irfftn(psi, axes = (0,1,2)).astype("float32")
    #CIC interpolation                                                                                                                                                                                                                                                                                                        


    vx = np.empty(len(pos)).astype("float32")
    vy = np.empty(len(pos)).astype("float32")
    vz = np.empty(len(pos)).astype("float32")
    #CIC interpolation scheme implemented in pylians (C) Francisco Villaescusa-Navarro
    MASL.CIC_interp(psi[:,:,:,0], box_size[0], pos+(box_size/2).astype("float32"), vx)
    MASL.CIC_interp(psi[:,:,:,1], box_size[0], pos+(box_size/2).astype("float32"), vy)
    MASL.CIC_interp(psi[:,:,:,2], box_size[0], pos+(box_size/2).astype("float32"), vz)

    del psi
    gc.collect()
    
    pos += box_size/2


    pos[:,0] = (pos[:,0] + vx)%box_size[0]
    pos[:,1] = (pos[:,1] + vy)%box_size[1]
    pos[:,2] = (pos[:,2] + vz)%box_size[2]

    pos -= box_size/2
    del vx, vy, vz
    gc.collect()
    return(pos)

