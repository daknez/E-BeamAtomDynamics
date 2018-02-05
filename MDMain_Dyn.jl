# ---------------------------------------------------------------------------
# S   T   A   R   T 
# ---------------------------------------------------------------------------

Pkg.add("MAT")	
Pkg.add("Combinatorics")
@everywhere using MAT
@everywhere using Combinatorics

@everywhere include("Functions.jl");	# Functions
@everywhere include("Constants.jl");	# Constants

# program files:
step1 = "Thermalize.jl"
step2 = "ThermalizeoncomplexSubstrate.jl"
exper = "Experiment_Dyn.jl"
main = "MDMain_Dyn.jl"

HT = 300;			# electron energy in keV
Temperature = 300;	# System temperature in K
TempR = 500;			# elevated relaxation temperature in K
element1 = "Au";
potential1 = "QSC";
#density = 8908; # Ni in kg/m³
#density = 10490; # Ag in kg/m³
density = 19320; #Au in kg/m³
nPart = 350;

# filenames:
postfix = string(element1,"_20180118_",nPart,"atoms_",potential1,"_Dynamics_",HT," kV_",Temperature,"K_TempR",TempR,"K")
infilename = string("random",element1,nPart);
prefix = "randsub_Lewis"
outfilename = string(prefix,postfix);

# create folder for results
if isdir(outfilename) .== false
	mkdir(outfilename)
end

# substrate parameters:

# random carbon arrangement:
LS = 4E-9;		# lateral size
dS = 0.4E-9;		# thickness
densityS = 2000;	# density

#Ni-C:
# Huang.2003, Ryu2010: 
# original:
# sigS = 2.852E-10;
# epsS = 0.023049 .* elch;

#Ag-C (Akbarzadeh.2014):
#epsS = 0.0301*elch; 
#sigS = 3.006e-10;

#Au-C: (Lewis.2000)
sigS = 2.74e-10;
epsS = 0.022*elch;

# Arcidiacono.2005:
# epsS = 0.01273*elch;
# sigS = 2.9943e-10;

# copy input files:
cp(string(element1,HT,".mat"),string(outfilename,"\\",element1,HT,".mat"),remove_destination=true)	# copy scattering cross section table	
cp(step1,string(outfilename,"\\",step1),remove_destination=true)					# copy program file for 1st equilibration
cp(step2,string(outfilename,"\\",step2),remove_destination=true)					# copy program file for 2nd equilibration			
cp(exper,string(outfilename,"\\",exper),remove_destination=true)					# copy program file for experiment
cp(main,string(outfilename,"\\",main),remove_destination=true)						# copy this file
cp("Functions.jl",string(outfilename,"\\","Functions.jl"),remove_destination=true)					
cp("Constants.jl",string(outfilename,"\\","Constants.jl"),remove_destination=true)
				
juliafolder = pwd();
cd(outfilename)

(a1, eps1, mSC1, nSC1, C1, m1) = SetMatParameters(element1,potential1);		# Get material parameters for given material and potential

rcluster = ((3*nPart*m1)/(4*pi*density))^(1/3);	# estimated cluster radius assuming a spherical geometry with given density and atom number

# === Set initial geometry ===s

#create random cluster and relax with MC:
(coords1, L) = create_rnd_cluster(nPart,density,m1);
nDim = size(coords1,1);

# --------- read geometry data from file --------------

# structfile = matopen(string(infilename,".mat"));
# coords = read(structfile, "coords");
# coords1 = coords[:,:,end];
# L = read(structfile, "L");
# #subz = read(structfile, "subz");
# close(structfile)
# nPart = size(coords1,2);
# nDim = size(coords1,1);
# postfix = string(nPart," atoms",postfix);

# # convert A to m:
# L = L * 1E-10;
# coords1 = coords1.*1E-10#.*0.8;

# MC equilibration:
trialsperatom = 3000;
sampleFreq = 1000;
startTemp = 1000;
# use simplified Metropolis-Monte-Carlo calculation to approximate aquilibrium structure:
coords1 = MC_Cluster_relaxation(coords1,L,element1,potential1,trialsperatom,sampleFreq,startTemp,outfilename);

#set centre of mass of the particle to 0 0 0:
centreofmass = sum(coords1,2)./nPart;
coords1 = broadcast(-,coords1,centreofmass);
coords1[3,:] = coords1[3,:]+2E-9;

#------------------------------------------
# creating substrate consisting of single atoms:
coordssub = create_aCsubstrate(LS,dS,densityS);
substratefilename = string("Csubstr_L",LS,"_d",dS,"_dens",densityS)
outfilename = string(prefix,infilename,substratefilename);
nPartS = size(coordssub,2);
centreofmassS = sum(coordssub,2)./nPartS;
coordssub[1:3,:] = broadcast(-,coordssub[1:3,:],centreofmassS[1:3]);
subz=0;
writeXYZ(coordssub.*1E10,"C",substratefilename);
println(string("initial number of particles: ",nPart));
println(string("substrate particles: ",nPartS));
#--------------------------------------------------------------------------

sqrt2 = sqrt(2)/2;
rc = 2*a1; 				# cut off distance
maxd = rc^2;				# cut of distance squared to save computation time
rNN = a1 * sqrt2*1.15;	# define maximum distance for nearest neighbours calculation
rNN2 = rNN*rNN;

# ===========================
# P A R A M E T E R S 
# ===========================

println(string("Element1: ",element1));
println(string("initial number of particles: ",nPart));

t1=t2=t3=0;

# ========================================================================
# ================= M O L E C U L A R	D Y N A M I C S =================
# ========================================================================

# =========================================
# T H E R M A L I S A T I O N 
# =========================================

tic();

nSteps = 5000; 						# no of thermalization timesteps
dt = 10E-15; 						# starting time step
Temp = 0.4;							# thermalisatio temperature
nu = 2E12;							# thermostat parameter
sampleFreq = 10;						# sampling frequency
include(step1);						# perform thermalisation
t1=toc();

#writeXYZ(coords1.*1E10,element1,string(postfix,"_afterinittherm"));

# ==========================================================
# T H E R M A L I S A T I O N   O N   S U B S T R A T E
# ==========================================================
tic();

# parameters for thermalization on substrate: 
nSteps = 5000; 						# no of thermalization timesteps
dt = 16E-15;       					# Integration time step (Villareal)
Temp = Temperature;					# Temperature for thermalization
nu = 2E12;							# Thermostat parameter for thermalization
sampleFreq = 10;						# sampling frequency for thermalization

include(step2);

t2=toc();

#writeXYZ(coords1.*1E10,element1,string(postfix,"_afterthermonsub"));
writeXYZ2el(coords1.*1E10,coordssub.*1E10,element1,"C",string(element1,nPart,"C",nPartS,potential1,"_afterthermonsub"));


# =========================================
# E X P E R I M E N T
# =========================================
# for the experiment: 
Eel = HT * 1E3 * elch;					# energy of incoming electrons in J
electronFrequ = 5000;						# Frequency of incident electrons (maximum time between 2 electron impacts)
nelectrons = 500000;
nSteps = 20000;     						# maximum number of simulation time steps per electron (in integration steps)
if element1 == "Au"
	dtE = 20E-15;     					# Integration time after Villarreal.2014
elseif element1 == "Ag"
	dtE = 16E-15; 
elseif element1 == "Ni"
	dtE = 14E-15;
else
	dtE = 10e-15;
end
Temp = Temperature;   				# Simulation temperature in K						# Temperature during relaxation
nu = 2E12;      						# Thermostat parameter - frequency of collisions with the heat bath (0.1)
sampleFreq = 20;    					# Sampling frequency

include(exper);
#writeXYZ(coords1.*1E10,element1,string(postfix,"_afterexperiment"));


print_with_color(:green,"Simulation finished...");
println(t1)
println(t2)
println(t3)
cd(juliafolder)