import numpy as np

from scipy.integrate import quad
from scipy.interpolate import CubicSpline

from settings import Omega_m, Omega_l, zmin, zmax, z_interp

##################################
# r to z conversion (r in Mpc/h) # ##################################                                                                                                                                                                                                                                                                      
def integrand(z):
    zp1 = 1+z
    hz = np.sqrt(Omega_m*zp1*zp1*zp1+Omega_l+(1-Omega_m-Omega_l)*zp1*zp1)
    return(2997.92458/hz)

z = np.linspace(zmin,zmax,z_interp)
r = np.zeros_like(z)

for i,red in enumerate(z):
    r[i] = quad(integrand, 0, red)[0]

r_to_z = CubicSpline(r,z)
z_to_r = CubicSpline(z,r)
