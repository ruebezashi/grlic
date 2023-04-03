import numpy as np

from settings import nbins, deg, n_provided, bufferr, buffergrid

# returns the polynomial function values of r, of degree deg with fitting parameters p
def f(r, p):

    f = 0
    for i in range(deg+1):
        f += r**i * p[deg-i]
    return(f)


#estimates n(r) from the list of comoving distances r by creating a number counts N(r) histogram, and fitting a polynomial to the n(r) derived from the survey shape
def fitted(r, mumin):
    
    #create number counts histogram of N(r)
    rhist = np.histogram(r, bins = nbins)
    counts, edges = rhist[0], rhist[1]
    deltabins = (edges[1] - edges[0])
    bins = (edges + deltabins/2)[:-1]

    #convert N(r) to n(r) taking the survey volume into account
    n_list = counts/bins**2/2/np.pi/(1-mumin)/deltabins
    r_list = bins

    #polyfit to the n_list(r_list), weighting by the inverse of the Poisson error within each shell
    upperbins = edges[1:]
    lowerbins = edges[:-1]
    volumes = 2*np.pi*(1-mumin)/3*(upperbins**3-lowerbins**3)
    p = np.polyfit(r_list, n_list, deg, w = volumes/np.sqrt(counts))
    n_list = f(r_list, p)

    return(r_list, n_list)


#reads in tabulated n(r) provided by user
def tabulated(rmin, rmax):

    r_list, n_list = np.loadtxt(n_provided, usecols = [0,1], unpack = True)
    mask = (r_list >= rmin) & (r_list <= rmax)
    
    return(r_list[mask], n_list[mask])


#extends the lists of comoving distances r_list and corresponding comoving number density n_list with a buffer zone of length bufferr
#returns buffered r and n arrays, as well as near and far constant number density
def buffered(r_list, n_list):
    
    #lower end buffer zone
    r0_list = np.linspace(0,r_list[0],buffergrid)
    n0_list = np.ones(buffergrid)*n_list[0]

    #append buffer to lower end
    r_list = np.append(r0_list[:-1], r_list)
    n_list = np.append(n0_list[:-1], n_list)

    #upper end buffer zone
    r1_list = np.linspace(r_list[-1], r_list[-1] + bufferr, buffergrid)
    n1_list = np.ones(buffergrid)*n_list[-1]

    #append buffer to upper end
    r_list = np.append(r_list, r1_list[1:])
    n_list = np.append(n_list, n1_list[1:])
    return(r_list, n_list, n0_list[-1], n1_list[0])

