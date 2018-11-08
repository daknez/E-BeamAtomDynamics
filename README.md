## E-BeamAtomDynamics

Collection of scripts for the simulation of electron beam induced dynamics in supported metallic clusters. The algorithm uses molecular dynamics, based on the Velocity-Verlet scheme to calculate trajectories of all atoms in a metallic cluster after elastic scattering with a swift electron. Scattering events are sampled randomly from NIST cross section tables. Data for Au and Ag for 300 keV electrons is stored in the corresponding Matlab data files. Cross section data for other elements and electron energies can be found at https://srdata.nist.gov/srd64/.

See our paper for details on the algorithm:

Knez, D. et al. Modelling electron beam induced dynamics in metallic nanoclusters. Ultramicroscopy 192, 69–79 (2018). DOI: 10.1016/j.ultramic.2018.05.007 (https://www.sciencedirect.com/science/article/pii/S0304399118300184)

Following files are needed to perform the simulation:
- MDMain.jl Generates the geometry, sets parameters and calls the single simulation steps.
- Thermalize.jl Equilibrates the cluster geometry.
- Thermalizeoncomplexsubstrate.jl The cluster is landed on an amorphous carbon substrate comprised of single atoms with fixed positions.
- Experiment.jl Performs the actual electron scattering experiment.
- Function.jl Contains all functions needed for the calculations.
- Constants.jl  Contains some constants.
- Ag300.mat and Au300.mat NIST differential cross section tables
