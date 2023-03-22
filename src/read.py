import numpy as np

from settings import data

#loads redshift, cosine of polar angle and azimuthal angle from catalogue
def zmuphi():
    print("reading data file...")
    z, mu, phi = np.loadtxt(data, usecols = [0,1,2], unpack = True, dtype = "float32")
    print("done")
    return(z,mu,phi)


#reads redshifts of catalogue
def z():
    print("reading redshifts of data file...")
    z = np.loadtxt(data, usecols = 0, dtype = "float32")
    print("done")
    return(z)

