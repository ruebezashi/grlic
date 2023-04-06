## grlic: glass-like random catalogues for cosmological light cones

## Description

grlic is a python tool which produces glass-like random catalogues with a specified radially-dependent number density distribution. Such catalogues can be used for precise estimates of N-point correlation functions of tracers on the backward light cone. The Landy-Szalay estimator of the two-point correlation function using these glass-like catalogues is unbiased and has much lower variance than the estimator using Poisson-sampled random catalogues.

## Requirements

grlic needs the following python packages:

 - Pylians3: https://pylians3.readthedocs.io/en/master/construction.html
 - numpy
 - scipy
 - numba

## Installation

In the current version, grlic can simply be run by cloning the repository and running:

python main.py

from inside the src directory.

This will create a glass-like random catalogue according to the settings specified in settings.py.


## Usage

Currently, there are two ways of generating a glass with grlic:

    1. The user can provide a path to a data catalogue in the settings.py file (data = ""), which contains three columns: the redshift (z), cosine angle of the polar angle (µ) and the azimuthal angle (phi) of each object in the catalogue. grlic will then generate a glass that follows the same distribution as the provided data catalogue.
    2. The user can also provide a file that contains the number density n at each comoving distance r (n_provided = ""), here the first column should be the comoving distance r and the second one should be the number density in units of (h/Mpc)^3. grlic will then generate a glass that follows the provided number density distribution.

Per default, an example is provided in the settings, which reads a tabulated n(r) distribution from "n_of_r.dat" and generates a glass saved as "glass.dat" within a survey volume with minimal comoving distance rmin, maximum comoving distance rmax and cosine of the opening angle mumin. A central slice through the glass, with a thickness of 50Mpc/h, can be plotted with the "plot_glass.py" script.

## Support
Any issues and questions will be answered by the author of the code, Sebastian Schulz:

E-Mail: sebastian.schulz@uzh.ch

## Data repository
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7799509.svg)](https://doi.org/10.5281/zenodo.7799509)



## Roadmap



## Contributing


## License

## Project status

