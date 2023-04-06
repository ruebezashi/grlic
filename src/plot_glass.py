import numpy as np
import matplotlib.pyplot as plt
import cosmology
import units


z, mu, phi = np.loadtxt("glass.dat",unpack = True)

r = cosmology.z_to_r(z)

x,y,z = units.spherical_to_cartesian(r,mu,phi)

mask = (y > -25) & (y < 25)

fix, ax = plt.subplots(1,1,figsize =(10,10))
ax.set_aspect("equal")
ax.scatter(x[mask], z[mask], s = 0.1)
plt.show()
