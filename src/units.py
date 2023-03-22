import numpy as np

def spherical_to_cartesian(r,mu,phi):
    x = r * np.sin(np.arccos(mu)) * np.cos(phi)
    y = r * np.sin(np.arccos(mu)) * np.sin(phi)
    z = r * mu
    return x,y,z

def cartesian_to_spherical(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    mu = z / r
    phi = np.arctan2(y, x)
    return r, mu, phi
    
