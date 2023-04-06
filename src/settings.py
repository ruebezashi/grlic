
# cosmological parameters ###########################

#density parameter of matter today
Omega_m= 0.318984543
#density parameter of dark energy today
Omega_l= 0.680947669


# survey parameters #################################

#cosine of the half-opening angle of the survey
mumin = 0.3

#minimum comoving distance, only used if tabulated n(r) is provided, otherwise this quantity is derived from the list of observed redshifts in the file
rmin = 200

#maxmimum comoving distance, only used if tabulated n(r) is provided, otherwise this quantity is derived from the list of observed redshifts in the file
rmax = 1300


# input files #######################################

#path to data catalogue, expects three columns: redshift z, angular position mu and phi, where mu is cos(theta), theta the polar angle and phi the azimuthal angle
data = "" 

#path to tabulated number density file, reads out first two columns in the file, r in Mpc/h and n in h^3 Mpc^(-3)
n_provided = "n_of_r.dat"

# output files ######################################

#desired path to the created glass catalogue
glass_name = "glass.dat"

# box grid parameters ###############################

#box mode, either cubic box with observer in center (useful for full-sky surveys at low redshifts), or cuboid box with observer outside of the box (for deep pencil beams)
#only available option is "spherical" ("pencil" to be implemented)
box_mode = "spherical"

#size of buffer zone in Mpc/h
bufferr = 400

#grid cell side length in Mpc/h
#will be adapted slightly to fit even number of cells into the box
cell_size = 10


# glass parameters #################################

#alpha factor that determines the number of objects in the random catalogue
alpha = 1

#number of Zeldovich iterations
n_iter = 2


#####################################################
# precision settings                                #
#####################################################

#minimum redshift of the cubic spline interpolation for translating between z and r
zmin = 2e-3
#maximum redshift of the cubic spline interpolation for translating between z and r
zmax = 3.998
#number of points in the z grid which is used in the z to r interpolation
z_interp = 1000

#number of bins in the histogram of n(r) used to perform the polynomial fit
nbins = 100
#degree of the polynomial fit to the number density histogram, default = 3, degrees from 1 to 6 are supported
deg = 3

#number of points in the buffer zone that is appended to the list of r values and comoving number densities inside the survey volume
buffergrid = 100

