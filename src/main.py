import matplotlib.pyplot as plt
# numpy library
import numpy as np

# garbage collection
import gc

# scipy library
from scipy.interpolate import CubicSpline
from scipy.integrate import quad

from settings import *
import read
import cosmology
import n_of_r
import box_particles
import poisson
import units
import zeldovich

def create_glass():
    if n_provided == "":
        #read file
        z = read.z()
        N = len(z)
        print("number of particles in the data catalogue: ", N)

        min_mu = mumin
        
        r = cosmology.z_to_r(z)
        min_r = np.min(r)
        max_r = np.max(r)
        print("minimum r: ", min_r, ", maximum r: ", max_r, ", minimum µ:", min_mu)
        r_list, n_list = n_of_r.fitted(r, min_mu)

    elif n_provided != "":
        print("using tabulated n(r)...")
        
        min_mu = mumin

        min_r = rmin
        max_r = rmax
        print("minimum r: ", min_r, ", maximum r: ", max_r, ", minimum µ:", min_mu)
        r_list, n_list = n_of_r.tabulated(min_r, max_r)

    r_list, n_list, n0, n1 = n_of_r.buffered(r_list,n_list)

    #create cubic spline interpolation for n(r)
    cs_n = CubicSpline(r_list,n_list)

    if n_provided != "":
        #integrate number density to get N
        N = int(quad(poisson.num_dens, rmin, rmax, args = (cs_n,  min_mu))[0])
        print("number of particles in the data catalogue, estimated from n(r): ", N)

    if box_mode == "spherical":

        #create box that encompasses full sphere plus buffer zone
        box_size = np.ones(3)*(max_r + bufferr)*2
        print("size of box: ",box_size)
        n_grid, r_grid = box_particles.n_grid_cube(box_size, cs_n)

        print("generating poisson random in survey...")
        #generate Poisson-sampled n(r) random within survey volume
        x_r, y_r, z_r = poisson.random_full(cs_n, N, min_r, max_r, min_mu)
        x_r, y_r, z_r = units.spherical_to_cartesian(x_r, y_r, z_r)

        print("generate poisson random outside of survey...")
        #generate Poisson-sampled uniform random outside of survey volume
        out_x, out_y, out_z = poisson.random_full_complement(min_r, n0, max_r, n1, box_size)
        
        pos = np.c_[np.append(x_r,out_x), np.append(y_r, out_y), np.append(z_r, out_z)].astype("float32")

        del x_r, y_r, z_r, out_x, out_y, out_z
        gc.collect()
        print("calculate density contrast")
        #calculate density contrast with respect to n(r) defined on n_grid
        delta = box_particles.delta_grid_cube(pos, box_size, n_grid)

        if n_iter == "":
            print("Performing Zeldovich iterations until the initial Poisson variance has been reduced by " + str(delta_reduction) + " percent ... \n")
            #calculate expected Poisson variance of CIC deltas in each grid cell
            ncells = (box_size/cell_size).astype("int") + (box_size/cell_size).astype("int")%2
            act_cell_size = box_size[0]/ncells[0]
            exp_N_per_cell = (n_grid*act_cell_size**3)
            exp_var = np.mean(exp_N_per_cell)*8/27
            act_var = np.mean((delta*exp_N_per_cell)**2)
            perc_reduction = (1 - act_var/exp_var)*100
            #iterate Zeldovich displacement
            iter_steps = 0
            while (perc_reduction < delta_reduction) & (iter_steps < n_iter_max):
                pos = zeldovich.displace_reverse(pos, delta, box_size)
                delta = box_particles.delta_grid_cube(pos, box_size, n_grid)
        
                perc_reduction = (1-np.mean((delta*exp_N_per_cell)**2)/exp_var)*100
                print("Poisson variance reduced by " + str(perc_reduction) + " percent \n")

                iter_steps += 1
                if iter_steps >= n_iter_max:
                    print("Maximum number of iterations reached, stopping iteration of Zeldovich displacement.")
        if n_iter != "":
            print("Performing " + str(n_iter) + " Zeldovich iterations... \n")
            #iterate Zeldovich displacement
            for i in range(n_iter):
                pos = zeldovich.displace_reverse(pos, delta, box_size)
                delta = box_particles.delta_grid_cube(pos, box_size, n_grid)
        
        #return r, mu, phi of random catalogue within survey volume
        pos[:,0], pos[:,1], pos[:,2] = units.cartesian_to_spherical(pos[:,0], pos[:,1], pos[:,2])
        
        mask = (pos[:,0] >= min_r) & (pos[:,0] <= max_r) & (pos[:,1] >= min_mu)
        pos[:,0][mask] = cosmology.r_to_z(pos[:,0][mask])
        
        np.savetxt(glass_name, pos[mask])

    '''
    elif box_mode == "pencil":
        #make survey volume a bit larger to keep edge effects away from survey edges
        min_r_buff = max(0, min_r-bufferr)
        max_r_buff = max_r + bufferr
        min_mu_buff = max(-1, min_mu - 0.1*min_mu)
        rnd_r, rnd_mu, rnd_phi = poisson.random_pencil(r_list, n_list, N, min_r_buff, max_r_buff, min_mu_buff)
    '''

create_glass()
