import numpy as np
import gc

import MAS_library as MASL

from settings import cell_size, alpha, mumin, bufferr


#assigns background number density n(r) to each grid cell in a box with dimensions box_size[0], box_size[1], box_size[2]
def n_grid_cube(box_size, n):

    #determine number of cells from cell_size and box_size, needs to be even for FFT
    ncells = (box_size/cell_size).astype("int") + (box_size/cell_size).astype("int")%2
    print("number of cells in grid: ", ncells, ", actual cell size: ", box_size/ncells)
    xi = ((np.arange(0,ncells[0])+0.5)*box_size[0]/ncells[0]-box_size[0]/2).astype("float32")
    yi = ((np.arange(0,ncells[1])+0.5)*box_size[1]/ncells[1]-box_size[1]/2).astype("float32")
    zi = ((np.arange(0,ncells[2])+0.5)*box_size[2]/ncells[2]-box_size[2]/2).astype("float32")
    print("calculate r for every grid cell")
    #calculate r for every grid cell 
    XI,YI,ZI = np.meshgrid(xi,yi,zi, indexing = "ij")
    del xi,yi,zi
    gc.collect()
    r = np.sqrt(XI**2 + YI**2 + ZI**2)
    del XI,YI,ZI
    gc.collect()
    n_grid = n(r)*alpha
    return(n_grid, r)


#computes delta in each grid cell from particle positions, using CIC scheme implemented in pylians (C) Francisco Villaescusa-Navarro
def delta_grid_cube(pos, box_size, n_grid):

    ncells = (box_size/cell_size).astype("int") + (box_size/cell_size).astype("int")%2

    # define 3D density field                                                                                                                                                                                                                                                                                                 
    delta = np.zeros((ncells[0],ncells[1],ncells[2]), dtype=np.float32)
  

    # construct 3D density field                                                                                                                                                                                                                                                                                              
    MASL.MA(pos+(box_size/2-box_size/2/ncells).astype("float32"), delta, box_size[0], MAS = "CIC")
    #average background number of particles in each grid cell                                                                                                                                                                                                                                                                 
    N_grid = n_grid * (box_size[0]/ncells[0])**3

    # at this point, delta contains the effective number of particles in each voxel                                                                                                                                                                                                                                           
    # now compute overdensity and density constrast                                                                                                                                                                                                                                                                           
    delta = delta/N_grid - 1

    return(delta)

#computes delta in each grid cell from particle positions, using CIC scheme implemented in pylians (C) Francisco Villaescusa-Navarro
def delta_grid_cube_test(pos, box_size, n_grid):

    ncells = (box_size/cell_size).astype("int") + (box_size/cell_size).astype("int")%2

    # define 3D density field                                                                                                                                                                                                                                                                                                 
    delta = np.zeros((ncells[0],ncells[1],ncells[2]), dtype=np.float32)
  

    # construct 3D density field                                                                                                                                                                                                                                                                                              
    MASL.MA(pos+(box_size/2-box_size/2/ncells).astype("float32"), delta, box_size[0], MAS = "NGP")
    #average background number of particles in each grid cell                                                                                                                                                                                                                                                                 
    N_grid = n_grid * (box_size[0]/ncells[0])**3

    # at this point, delta contains the effective number of particles in each voxel                                                                                                                                                                                                                                           
    # now compute overdensity and density constrast                                                                                                                                                                                                                                                                           
    delta = delta/N_grid - 1

    return(delta)

'''
#assigns background number density n(r) to each grid cell in a box with dimensions box_size[0], box_size[1], box_size[2]
def n_grid_cuboid(box_size, n, obs):

    #determine number of cells from cell_size and box_size, needs to be even for FFT
    ncells = int(box_size/cell_size) + int(box_size/cell_size)%2
    xi = ((np.arange(0,ncells[0])+0.5)*box_size[0]/ncells[0]-obs[0]).astype("float32")
    yi = ((np.arange(0,ncells[1])+0.5)*box_size[1]/ncells[1]-obs[1]).astype("float32")
    zi = ((np.arange(0,ncells[2])+0.5)*box_size[2]/ncells[2]-obs[2]).astype("float32")
    #calculate r for every grid cell 
    XI,YI,ZI = np.meshgrid(xi,yi,zi, indexing = "ij")
    del xi,yi,zi
    gc.collect()
    r = np.sqrt(XI**2 + YI**2 + ZI**2)
    #calculte mu for every grid cell
    MU = Z/R
    MU[MU < mumin] = 0
    del XI,YI,ZI
    gc.collect()
    n_grid = n(r)*alpha
    return(n_grid)
'''
