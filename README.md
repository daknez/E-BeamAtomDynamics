## E-BeamAtomDynamics

Collection of scripts for the simulation of electron beam induced dynamics in supported metallic clusters. The algorithm uses molecular dynamics, based on the Velocity-Verlet scheme to calculate trajectories of all atoms in a metallic cluster after elastic scattering with a swift electron. Scattering events are sampled randomly from NIST cross section tables. Data for Au and Ag for 300 keV electrons is stored in the corresponding Matlab data files. Cross section data for other elements and electron energies can be found at https://srdata.nist.gov/srd64/.

See our paper for details on the algorithm:

Knez, D. et al. Modelling electron beam induced dynamics in metallic nanoclusters. Ultramicroscopy 192, 69–79 (2018). DOI: 10.1016/j.ultramic.2018.05.007 (https://www.sciencedirect.com/science/article/pii/S0304399118300184)

The repository contains following files:
- MDMain.jl Generates the geometry, sets parameters and calls the single simulation steps.
- 
