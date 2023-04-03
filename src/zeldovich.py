import gc
import numpy as np
import warnings

from scipy.fft import fftfreq
from scipy.fft import fft
from scipy.fft import ifft

import MAS_library as MASL

from settings import cell_size

def displace_reverse(pos, delta, box_size):


    ncells = (box_size/cell_size).astype("int") + (box_size/cell_size).astype("int")%2

    # fourier transform delta in 3D                                                                                                                                                                                                                                                                                           
    deltaf = fft(delta, axis = 0).astype("csingle")
    del delta
    deltaf = fft(deltaf, axis = 1).astype("csingle")
    deltaf = fft(deltaf, axis = 2).astype("csingle")

    #unit of fftfreq is in cycles / Mpc/h -> to get to h/Mpc multiply by 2pi                                                                                                                                                                                                                                                  
    kx = (fftfreq(ncells[0], box_size[0]/ncells[0])*2*np.pi).astype("csingle")
    ky = (fftfreq(ncells[1], box_size[1]/ncells[1])*2*np.pi).astype("csingle")
    kz = (fftfreq(ncells[2], box_size[2]/ncells[2])*2*np.pi).astype("csingle")
    #get FT-displacement field from FT-delta                                                                                                                                                                                                                                                                                  
    psi = np.zeros((ncells[0],ncells[1],ncells[2],3)).astype("csingle")
    KX, KY, KZ = np.meshgrid(kx,ky,kz, indexing = "ij")
    KSQ = (KX**2 + KY**2 + KZ**2)
    KSQ_st = np.stack((KSQ,KSQ,KSQ), axis = 3)
    KAR = np.stack((KX,KY,KZ),axis = 3)
    del KX, KY, KZ, KSQ, kx, ky, kz
    gc.collect()
    deltaf_st = np.stack((deltaf, deltaf, deltaf), axis = 3)
    del deltaf
    gc.collect()
    
    #reverse displacement field
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        psi = -1j*KAR / KSQ_st * deltaf_st
    psi[0,0,0] = 0
    psi[int(ncells[0]/2),:,:] = 0
    psi[:,int(ncells[1]/2),:] = 0
    psi[:,:,int(ncells[2]/2)] = 0

    #inverse FT of displacement Field                                                                                                                                                                                                                                                                                         
    ipsif = ifft(psi, axis = 0)
    del psi
    gc.collect()
    ipsif = ifft(ipsif, axis = 1)
    ipsif = ifft(ipsif, axis = 2)

    #CIC interpolation                                                                                                                                                                                                                                                                                                        
    ipsifx = ipsif[:,:,:,0].real
    ipsify = ipsif[:,:,:,1].real
    ipsifz = ipsif[:,:,:,2].real

    del ipsif
    gc.collect()

    vx = np.zeros(len(pos)).astype("float32")
    vy = np.zeros(len(pos)).astype("float32")
    vz = np.zeros(len(pos)).astype("float32")

    #CIC interpolation scheme implemented in pylians (C) Francisco Villaescusa-Navarro
    MASL.CIC_interp(ipsifx, box_size[0], pos+(box_size/2).astype("float32"), vx)
    MASL.CIC_interp(ipsify, box_size[0], pos+(box_size/2).astype("float32"), vy)
    MASL.CIC_interp(ipsifz, box_size[0], pos+(box_size/2).astype("float32"), vz)

    del ipsifx
    del ipsify
    del ipsifz
    gc.collect()
    
    pos += box_size/2


    pos[:,0] = (pos[:,0] + vx)%box_size[0]
    pos[:,1] = (pos[:,1] + vy)%box_size[1]
    pos[:,2] = (pos[:,2] + vz)%box_size[2]

    pos -= box_size/2
    del vx, vy, vz
    gc.collect()
    return(pos)

