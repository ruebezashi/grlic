
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad

from settings import alpha, deg

#given the n(r), returns the probability density taking the survey volume into account (normalizing by the number of objects in the data catalogue N)
def prob_dens(r, n, N, mumin):
    
    return(2*(1-mumin)*np.pi/N*n(r)*r**2)


#integrand over which the integral to get the number of objects within the survey N is performed
def num_dens(r, n, mumin):
    
    return(2*(1-mumin)*np.pi*n(r)*r**2)
'''
#returns the cumulative probability given a polynomial with coefficients p, taking the survey volume into account (normalizing by the number of objects in the data catalogue N)
def cumul_polynomial(r, p, N, mumin, rmin):
    
    f = 0
    for i in range(deg+1):
        f += 2*(1-mumin)*np.pi/N *((r**(i+3) * p[deg-i] / (i+3)) - (rmin**(i+3) * p[deg-i] / (i+3)))
    return(f)
'''

#generates a poisson sampled random catalogue on the full spherical shell with alpha times N_cat objects in the actual survey volume, based on the interpolated n(r) (using the cs function) of a data catalogue with N_cat objects, assuming a spherical cone shape with half-opening cosine angle mu_min and limits r_min to r_max
#the survey volume will contain alpha times N_cat objects, the rest of the spherical shell will be filled up with objects according to the same n(r)
def random_full(cs, N_cat, rmin, rmax, mumin):
    
    #10.000 interpolated points
    x = np.linspace(rmin, rmax, 10000)
    
    #initiate cumulative probability array
    cumul = np.zeros_like(x)

    #integrate probability density
    for i in range(len(x)):
        cumul[i] = quad(prob_dens, rmin, x[i], args = (cs, N_cat, mumin))[0]
    
    #invert cumulative probability with cubic spline interpolation
    inverse_cumul = CubicSpline(cumul, x)
    
    #sample from inverse cumulative probability randomly
    sample = np.random.rand(alpha * N_cat)
    #convert sample to distances
    rand_r = inverse_cumul(sample)
    
    #generate random mus
    rand_mu = np.random.rand(alpha * N_cat)*(1-mumin)+mumin

    #generate random phis
    rand_phi = np.random.rand(alpha * N_cat)*(2*np.pi)

    # fill up the rest of the spherical shell volume with same n(r)
    ratio = (1+mumin)/(1-mumin)
    N_noncat = int( ratio * N_cat )
    rand_r2 = inverse_cumul(np.random.rand(alpha * N_noncat))
    rand_mu2 = np.random.rand(alpha * N_noncat) * (-1 - mumin) + mumin
    rand_phi2 = np.random.rand(alpha * N_noncat) * 2 * np.pi

    rand_r = np.append(rand_r, rand_r2)
    rand_mu = np.append(rand_mu, rand_mu2)
    rand_phi = np.append(rand_phi, rand_phi2)
    
    return rand_r, rand_mu, rand_phi



#complements the poisson sampled spherical shell random catalogue in the box: below rmin we reproduce nmin and above rmax we reproduce nmax
def random_full_complement(rmin, n0, rmax, n1, box_size):
    #inner region
    x_i = np.random.rand(int(n0*box_size[0]**3*alpha)).astype("float32") * box_size[0] - box_size[0]/2
    y_i = np.random.rand(int(n0*box_size[1]**3*alpha)).astype("float32") * box_size[1] - box_size[1]/2
    z_i = np.random.rand(int(n0*box_size[2]**3*alpha)).astype("float32") * box_size[2] - box_size[2]/2
    r = np.sqrt(x_i**2 + y_i**2 + z_i**2)
    x_i = x_i[r < rmin]
    y_i = y_i[r < rmin]
    z_i = z_i[r < rmin]
    
    #outer region
    x_o = np.random.rand(int(n1*box_size[0]**3*alpha)).astype("float32") * box_size[0] - box_size[0]/2
    y_o = np.random.rand(int(n1*box_size[1]**3*alpha)).astype("float32") * box_size[1] - box_size[1]/2
    z_o = np.random.rand(int(n1*box_size[2]**3*alpha)).astype("float32") * box_size[2] - box_size[2]/2
    r = np.sqrt(x_o**2 + y_o**2 + z_o**2)
    x_o = x_o[r > rmax]
    y_o = y_o[r > rmax]
    z_o = z_o[r > rmax]
    return np.append(x_i,x_o), np.append(y_i,y_o), np.append(z_i,z_o)
    
'''
#generates a poisson sampled random catalogue with alpha times N_cat objects, based on the n(r) of a data catalogue with N_cat objects, assuming a spherical cone shape with half-opening cosine angle mu_min and limits r_min to r_max
def random_pencil(r, n, N_cat, rmin, rmax, mumin):


    #cubic spline interpolate n(r) finely for accurate numerical integration
    spline1 = CubicSpline(r,n)
    
    #10.000 interpolated points
    x = np.linspace(rmin, rmax, 10000)
    
    #initiate cumulative probability array
    cumul = np.zeros_like(x)
    
    #integrate probability density
    for i in range(len(x)):
        cumul[i] = quad(prob_dens, rmin, x[i], args = (spline1, N_cat, mumin))[0]
    
    #invert cumulative probability with cubic spline interpolation
    inverse_cumul = CubicSpline(cumul, x)
    
    #sample from inverse cumulative probability randomly
    sample = np.random.rand(alpha * N_cat)
    #convert sample to distances
    rand_r = inverse_cumul(sample)
    
    #generate random mus
    rand_mu = np.random.rand(alpha * N_cat)*(1-mumin)+mumin
    #generate random phis
    rand_phi = np.random.rand(alpha * N_cat)*(2*np.pi)
    
    return rand_r, rand_mu, rand_phi 
'''






