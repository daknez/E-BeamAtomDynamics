	
	scatterfile = matopen(string(element1,HT,".mat"));	# open Matlab file containing the tabulated scattering cross section values
	scattercross = read(scatterfile, "cross_section");	# list of cross sections in units of a0^2
	scatterangle = read(scatterfile, "deg");			# list of angles in deg
	close(scatterfile)
	println(string(scatterfile," loaded."));
	
	dt = dtE;
	dt2 = dt*dt;
	timeC = dt*maximum((forces1[1,:].^2+forces1[2,:].^2+forces1[3,:].^2).^0.5);	# after Marks.2005
	
	time = 0; # Following simulation time
	
	samplecount = 1;
	sputterindex = 1;
	
	displThr = a1.*(2^0.5/2)*0.98;	# displacement threshold for detection of a displaced atom
		
	nPartStart = nPart;
	
	maxnrofpairs = size(collect(combinations(1:nPart,2)),1);
	
	#experimental parameters:
	Et = 0.0;
	Etmax = 2*Eel*(Eel+2*Eel0)/(c^2*m1);
	Etmin = 0.4*elch;
	nel = 1;
	partdisp = 0.0;
	initcoords = coords1;
	precoords = zeros(size(coords1,2));
	indextemp = 0;
	
	# set theta range and power law-parameters (300 kV):
	#thetamin = 0;
	thetamin = 2*asin(sqrt(Etmin/Etmax))*180/pi;
	thetamax = 180;
	
	X = scatterangle[scatterangle.>=thetamin];
	X = X[X.<=thetamax];
	indmin = find(scatterangle.==minimum(X));
	indmax = find(scatterangle.==maximum(X));
	P = scattercross[indmin[1]:indmax[1]];

	#crosssec = 2*pi*sum(P.*sin(scatterangle[indmin[1]:indmax[1]]./180.*pi));
	crosstot = 2*pi*trapz(scatterangle./180.*pi,scattercross.*sin(scatterangle./180.*pi))*a0^2;	# trapz performs a numerical integration after the trapezoidal rule
	crosssec = 2*pi*trapz(X./180.*pi,P.*sin(scatterangle[indmin[1]:indmax[1]]./180.*pi))*a0^2;
	
	scalingfactor = crosstot/crosssec;
	arealdensity = 350/5e-18;	# estimated projected area for a cluster comprised of 350 atoms: 5 nm²
	
	shottimestep = 50;
	waitsteps = 250;				# steps between scattering event and 1st displacement check (500)
	secdispcheckstep = 50;		# steps between 1st positive displacement check and 2nd displacement check (100)
	maxelFreq = electronFrequ;
	minelFreq = 700;
	
	shottime = 0.2e-12; #shottimestep*dt;
	waittime = 3e-12 #waitsteps*dt;
	secdispchecktime = 1e-12 # secdispcheckstep*dt;
	equilibrationtime = 50e-12 # electronFrequ*dt;
	
	println("-----------------------------");
	print_with_color(:green,"Experiment started\n");
	println(string("Temp=",Temp,"K, nu=",nu/1E12," THz", ", Time step: ",round(dt*1E15,1)," fs"));
	println("-----------------------------");

	shotelectrons_temp = zeros(7,nelectrons);
	sputteringevents = zeros(6,nelectrons);
	nel=1;
	theta = 0;
	atomnr = 0;
	resetconfig = false;			# reset configuration every electron?
	fulllog = false;			# write full log at every displacement event (CAUTION: Creates a lot of data!)
	
	index = 0;
	
	while nel < nelectrons && nPart >= 5
		
		tic();
		samplecount = 1;
		displacement = false;
		dispcheck1 = false;
		thermostaton = true;
		println("Thermostat on")
		shotflag = false;
		dispcheckflag = false;

		timedispevent = 0;

		nPart = size(coords1,2);
		
		countNN1 = zeros(nPart);

		# calculate initial thermal velocities:
		velSTDT1 = sqrt(kB*Temp/m1);
		vels1 = InitVels(nPart, nDim, velSTDT1, Temp);
		
		# prepare variables for sampling:
		coord_temp = zeros(nDim,nPartStart,10*convert(Int64,round((maxelFreq+shottimestep)/sampleFreq)));
		vels_temp = zeros(nDim,nPartStart,10*convert(Int64,round((maxelFreq+shottimestep)/sampleFreq)));
		countNN_temp = zeros(nPartStart,10*convert(Int64,round(nSteps/sampleFreq)+1));
		time_temp = zeros(1,10*convert(Int64,round((maxelFreq+shottimestep)/sampleFreq)));
		E_temp = zeros(2,10*convert(Int64,round((maxelFreq+shottimestep)/sampleFreq)));
		
		# Calculate initial forces
		forces1 = zeros(size(coords1));
		
		(rhos,pairs) = calculate_rhopairs(coords1,a1, mSC1, maxd, maxnrofpairs)
		forces1 = SC_ForcePL(coords1,rhos,L,pairs,a1,eps1,mSC1,nSC1,C1, maxd);

		print_with_color(:green,"next electron nr: ");
		println(nel);
		print_with_color(:green,"number of electrons in reality:");
		print_with_color(:green,"Nr of atoms: ");
		println(nPart);
		thermalized = false;
	
		
		for stept=1:nSteps

			tempcount = 0;
			# ---------------------------------
			# shoot an electron:
			# ---------------------------------
			if (timedispevent >= shottime) && shotflag == false
				# Output variables:
				(Ekin1, Epot1, countNN1) = calcEnergies(coords1,L,maxd,a1,eps1,m1,mSC1,nSC1,C1,vels1,rNN2);
				
				Ekin_sum = sum(Ekin1)
				Epot_sum = sum(Epot1)

				println(string("step no: ",stept,"; number of particles: ",nPart))
				println(string("Temperature of Andersen Thermostat: ",Temp," K"))
				println(string("number of atoms considered for thermostat: ",tempcount))
				println(string("elapsed time in simulation: ",time," s"))
				println(string("elapsed time in current displacement event: ",timedispevent," s"))
				println(string("current time step: ",dt," s"))
				println(string("Ekin = ",Ekin_sum," J = ",Ekin_sum/elch," eV"))

				# Temperatur in non-periodic system 2*Ekin/3N-6 (3f translation & 3 f rotation)
				println(string("Temperature from Ekin = ",2.*Ekin_sum/(3.*nPart-6)/kB," K"))
				println(string("Epot = ",Epot_sum," J = ",Epot_sum/elch," eV"))
				println(string("Thermostat is on: ", thermostaton))
				println("-------------------------------------------------------------")
				
				precoords = zeros(size(coords1,2));
				precoords = coords1; # save atom coordinates displacement check
				
				# ------------- choose theta (scattering angle) ------------------
				
				#theta from NIST-file:
				cdf = cumsum(P);
				theta = X[sum(cdf.< rand()*cdf[end]) + 1];			# Generating a non uniform discrete random variable from, scattering data
				println(string("Scattering angle chosen: ",theta))
				theta = theta*pi/180;								# Convert scattering angle in rad

				#theta equally distributed:
				#theta = rand()*thetamax+thetamin;
				#theta = theta*pi/180;
				#-----------------------------------------------------------------
				atomnr = rand(1:nPart);
				(v_sc, Et, psi, phi) = calcscatteringvelvector(m1,Eel,theta);
				
				# calculate new velocity vector of atomnr:
				#v_sc = [v_sc[1],v_sc[3],v_sc[2]]			# beam in y direction instead of z (electron beam from the side)
				vels1[:,atomnr] = vels1[:,atomnr] -  v_sc	# supported particle configuration
				#vels1[:,atomnr] = vels1[:,atomnr] +  v_sc 	# free hanging particle configuration
	
				# Save parameters to variable:
				shotelectrons_temp[6,nel] = countNN1[atomnr];
				shotelectrons_temp[7,nel] = 1;
				
				shotelectrons_temp[1,nel] = atomnr;
				shotelectrons_temp[2,nel] = Et;
				shotelectrons_temp[5,nel] = theta;
				
				thermostaton = false;			# deactivate thermostat to avoid pertubations during scattering event
				println("Thermostat turned off")
				velSTDT1 = sqrt(kB*TempR/m1);  # use different (higher) temperature TempR for relaxation when thermostat is reactivated
				
				shotflag = true;
				
				print_with_color(:green,"electron shot...\nEnergy transferred by electron [eV]: ");
				println(Et/elch)
				print_with_color(:green,"to atom number: ");
				println(atomnr)
			
			# ---------------------------------
			# Check for displacements (1st):
			# ---------------------------------

			elseif (timedispevent >= shottime+waittime) && shotflag && dispcheck1 == false && dispcheckflag == false
				
				dcoords = precoords - coords1;

				dcoordsabs = (dcoords[1,:].*dcoords[1,:]+dcoords[2,:].*dcoords[2,:]+dcoords[3,:].*dcoords[3,:]).^0.5;
				dsorted = sort(dcoordsabs);
				
				noofdisplacements = size(find((dcoordsabs .> displThr) .== true),1)	# number of atoms displaced over threshold
				index = find(dcoordsabs .== dsorted[end])		# find largest value
				index2 = find(dcoordsabs .== dsorted[end-1])	# find second largest value
				
				print_with_color(:green,"check atom displacements...\n maximum displacement found: ");
				println(string(dsorted[end]*1E9," nm (displacement threshold: ",displThr*1E9," nm) second displacement: ",dsorted[end-1]*1E9," nm"));
				
				print_with_color(:green,"from atom number: ");
				println(index[1]);
				
				# ---------------------------------
				# 1. check for atom displacements:
				# ---------------------------------
				if noofdisplacements >= 1   # atom movement of about 0.9 .* sigma is considered as permanent displacement due to electron impact
					
					println(string("calculated number of steps for current electron: ",stept));
					dispcheck1 = true;
					print_with_color(:yellow,string("displaced atom detected (",noofdisplacements," found)"));

				else
					coords1 = precoords;	# go back to initial coordinates
					dispcheckflag == true
					print_with_color(:red,"atom not displaced --> new electron!\n");
					break;
				end
				
			# ---------------------------------
			# Check for displacements (2nd):
			# ---------------------------------
			
			elseif dispcheck1 && displacement == false && timedispevent > (shottime + waittime + secdispchecktime)
				# when a displacement is detected in the first check, a second check some time after the first one is done
				print_with_color(:yellow,", check again...\n");
				dcoords = precoords - coords1;

				dcoordsabs = (dcoords[1,:].*dcoords[1,:]+dcoords[2,:].*dcoords[2,:]+dcoords[3,:].*dcoords[3,:]).^0.5;
				dsorted = sort(dcoordsabs);
				
				noofdisplacements = size(find((dcoordsabs .> displThr) .== true),1)	# number of atoms displaced over threshold
				indexnew = find(dcoordsabs .== dsorted[end])		# find largest value
				indexnew2 = find(dcoordsabs .== dsorted[end-1])	# find second largest value
				
				println(string("maximum displacement found: ",dsorted[end]*1E9," nm (displacement threshold: ",displThr*1E9," nm) second displacement found: ",dsorted[end-1]*1E9," nm"));
				print_with_color(:green,"from atom number: ");
				println(indexnew[1]);
				
				# ---------------------------------
				# 2. check for atom displacements:
				# ---------------------------------

				if size(find((dcoordsabs .> displThr) .== true),1) >= 1
				
					print_with_color(:red,string("atoms still displaced (", noofdisplacements," found)\n"));
					print_with_color(:red,"displacement confirmed.\n");
					println(string("calculated number of steps for current electron: ",stept));
					displacement = true;

					println("calculate relaxation:");
				else
					coords1 = precoords;	# go back to initial coordinates
					
					displacement = false;
					#dispcheck1 = false;
					print_with_color(:red,"displacement not confirmed. --> new electron!\n");
					break;
				end
				
				dispcheck1 = false
				thermostaton = true
				println("Thermostat on")
				dispcheckflag = true
				
				shotelectrons_temp[3,nel] = maximum(dcoordsabs);		
				shotelectrons_temp[4,nel] = index[1];
				
			elseif displacement && dispcheckflag && (timedispevent > (shottime + equilibrationtime))	# wait until relaxation and end
				
				print_with_color(:blue,"relaxation finished --> save results...\n");
				break;
			end


			# === 1st  integration step ===
			# Update positions:
			coords1 = coords1 + dt.*vels1 + 0.5.*dt2.*forces1./m1; 
				
			# Update velocities:
			vels1 = vels1 + 0.5.*dt.*forces1./m1;

			if mod(stept,20)==0
				if displacement == true
				
					# Check for lost atoms:
					#calculate center of particle (without accounting for mass):
					centreofmass = sum(coords1,2)./nPart;
					
					#-----------------------
					# ELEMENT 1:
					#-----------------------
					part=1;
					while part <= nPart
						# check if particle is outside of z-boundary
						if coords1[3,part] > L
							partdisp = part;
							lastcoords = coords1[:,partdisp];
							#delete particle outside of boundary:
							coords1 = coords1[:,[1:partdisp-1;partdisp+1:end]];
							vels1 = vels1[:,[1:partdisp-1;partdisp+1:end]];
							forces1 = forces1[:,[1:partdisp-1;partdisp+1:end]];
							countNN1 = countNN1[[1:partdisp-1;partdisp+1:end]];
							
							sputteringevents[1,sputterindex] = nel;		# electron number
							sputteringevents[2,sputterindex] = theta;	# scattering angle chosen
							sputteringevents[3,sputterindex] = Et;		# transferred Energy to atom
							sputteringevents[4,sputterindex] = atomnr;	# atom chosen for energy transfer
							sputteringevents[5,sputterindex] = partdisp;	# atom lost
							sputteringevents[6,sputterindex] = 1;		# element1
							nPart = size(coords1,2);	
							
							print_with_color(:red,string("one atom left the box (z-boundary, Last atom coordinates:",lastcoords,") atoms remaining: "));
							println(nPart)
							partdisp = 0;
							
							part = part - 1;
							sputterindex = sputterindex + 1;
						
						elseif (coords1[2,part]-centreofmass[2]) > L || (coords1[2,part]-centreofmass[2]) <-L
							partdisp = part;
							lastcoords = coords1[:,partdisp];
							#delete particle outside of boundary:
							coords1 = coords1[:,[1:partdisp-1;partdisp+1:end]];
							vels1 = vels1[:,[1:partdisp-1;partdisp+1:end]];
							forces1 = forces1[:,[1:partdisp-1;partdisp+1:end]];
							countNN1 = countNN1[[1:partdisp-1;partdisp+1:end]];
							
							sputteringevents[1,sputterindex] = nel;		# electron number
							sputteringevents[2,sputterindex] = theta;	# scattering angle chosen
							sputteringevents[3,sputterindex] = Et;		# transferred Energy to atom
							sputteringevents[4,sputterindex] = atomnr;	# atom chosen for energy transfer
							sputteringevents[5,sputterindex] = partdisp;	# atom lost
							sputteringevents[6,sputterindex] = 1;		# element1
							nPart = size(coords1,2);	
							
							print_with_color(:red,string("one atom left the box (y-boundary, Last atom coordinates:",lastcoords,") atoms remaining: "));
							println(nPart)
							partdisp = 0;
							
							part = part-1;
							sputterindex = sputterindex + 1;
							
						elseif (coords1[1,part]-centreofmass[1]) > L || (coords1[1,part]-centreofmass[1]) <-L
							partdisp = part;
							lastcoords = coords1[:,partdisp];
							#delete particle outside of boundary:
							coords1 = coords1[:,[1:partdisp-1;partdisp+1:end]];
							vels1 = vels1[:,[1:partdisp-1;partdisp+1:end]];
							forces1 = forces1[:,[1:partdisp-1;partdisp+1:end]];
							countNN1 = countNN1[[1:partdisp-1;partdisp+1:end]];
							
							sputteringevents[1,sputterindex] = nel;		# electron number
							sputteringevents[2,sputterindex] = theta;	# scattering angle chosen
							sputteringevents[3,sputterindex] = Et;		# transferred Energy to atom
							sputteringevents[4,sputterindex] = atomnr;	# atom chosen for energy transfer
							sputteringevents[5,sputterindex] = partdisp;	# atom lost
							sputteringevents[6,sputterindex] = 1;		# element1
							nPart = size(coords1,2);
							
							print_with_color(:red,string("one atom left the box (x-boundary, Last atom coordinates:",lastcoords,") atoms remaining: "));
							println(nPart);
							partdisp = 0;
							
							part = part-1;
							sputterindex = sputterindex + 1;
						end
						part = part + 1;
					end
				end
			end

		# Calculate new forces:
		
		# parallel:

		parad1 = calculate_rhopairs(coords1,a1, mSC1, maxd, maxnrofpairs);

		# substrate interaction
		parad2 = @spawn sub_LJforces(coordssub[:,1:convert(Int64,round(nPartS/3))],coords1,sigS,epsS);
		parad3 = @spawn sub_LJforces(coordssub[:,convert(Int64,round(nPartS/3))+1:convert(Int64,round(2*nPartS/3))],coords1,sigS,epsS);
		parad4 = @spawn sub_LJforces(coordssub[:,convert(Int64,round(2*nPartS/3))+1:end],coords1,sigS,epsS);
		
		(subforces1,thermpartlist1) = fetch(parad2);
		(subforces2,thermpartlist2) = fetch(parad3);
		(subforces3,thermpartlist3) = fetch(parad4);
		(rhos,pairs) = fetch(parad1);
				
		subforces = subforces1 + subforces2+subforces3;
		thermpartlist = thermpartlist1 | thermpartlist2 | thermpartlist3;	# list of atoms considered for surface thermostat

		forces1 = SC_ForcePL(coords1,rhos,L,pairs,a1,eps1,mSC1,nSC1,C1, maxd) + subforces;
			
			# Implement the surface Andersen thermostat
			if thermostaton
				list = find(thermpartlist.==true)
				for part = list
					# Test for collisions with the Andersen heat bath
					if rand() < nu.*dt
						vels1[:,part] = randGauss(0,velSTDT1,nDim);	
						vels1[3,part] = abs(vels1[3,part]);
					end
					tempcount = tempcount + 1;
				end
			end
			
			# === 2nd integration step ===
		  
			# Update velocities:
			vels1 = vels1 + 0.5.*dt.*forces1./m1;
			
			# Move time forward:
			time = time + dt;
			timedispevent = timedispevent + dt;
			
			# adaptive dt after Marks.2015:
			dt = timeC/maximum((forces1[1,:].^2+forces1[2,:].^2+forces1[3,:].^2).^0.5);
			dt2 = dt*dt;

			# === Sample ===
			if mod(stept,sampleFreq) == 0
				
			
				(Ekin1, Epot1, countNN1) = calcEnergies(coords1,L,maxd,a1,eps1,m1,mSC1,nSC1,C1,vels1,rNN2);
				
				Ekin_sum = sum(Ekin1)
				Epot_sum = sum(Epot1)
					
				if mod(samplecount,100) == 0	

					println(string("step no: ",stept,"; number of particles: ",nPart))
					println(string("Temperature of Andersen Thermostat: ",Temp," K"))
					println(string("number of atoms considered for thermostat: ",tempcount))
					println(string("elapsed time in simulation: ",time," s"))
					println(string("elapsed time in current displacement event: ",timedispevent," s"))
					println(string("current time step: ",dt," s"))
					println(string("Ekin = ",Ekin_sum," J = ",Ekin_sum/elch," eV"))
					# Temperatur in non-periodic system 2*Ekin/3N-6 (3f translation & 3 f rotation)
					println(string("Temperature from Ekin = ",2.*Ekin_sum/(3.*nPart-6)/kB," K"))
					println(string("Epot = ",Epot_sum," J = ",Epot_sum/elch," eV"))
					println(string("Thermostat is on: ", thermostaton))
					println("-------------------------------------------------------------")
				end
				
				# convert to SI units and save
				fillzero = zeros(3,size(coord_temp,2)-nPart);

				time_temp[1,samplecount] = time; 
				coord_temp[1:3,:,samplecount] = hcat(coords1,fillzero);
				vels_temp[1:3,:,samplecount] = hcat(vels1,fillzero);
				E_temp[1,samplecount] = Ekin_sum;
				E_temp[2,samplecount] = Epot_sum;
				
				countNN_temp[:,samplecount] = vcat(countNN1,zeros(size(coord_temp,2)-size(coords1,2)));
				samplecount = samplecount+1;
			end
		end

		if displacement
			# are there atoms with 0 coordination?
			countNNzero1 = find(countNN1.==0);
			if length(countNNzero1) > 0
				for ind = 1:length(countNNzero1)
					partdisp = countNNzero1[ind];
					coords1 = coords1[:,[1:partdisp-1;partdisp+1:end]];
					vels1 = vels1[:,[1:partdisp-1;partdisp+1:end]];
					forces1 = forces1[:,[1:partdisp-1;partdisp+1:end]];
					
					sputteringevents[1,sputterindex] = nel;		# electron number
					sputteringevents[2,sputterindex] = theta;	# scattering angle chosen
					sputteringevents[3,sputterindex] = Et;		# transferred Energy to atom
					sputteringevents[4,sputterindex] = atomnr;	# atom chosen for energy transfer
					sputteringevents[5,sputterindex] = partdisp;	# atom lost
					sputteringevents[6,sputterindex] = 1;		# element1
					
					sputterindex = sputterindex + 1;
				end
				nPart = size(coords1,2);				
					
				print_with_color(:red,string("one atom found without neighbours (atoms remaining from element1: "));
				println(string(nPart,")"))
				partdisp = 0;
			end
		end
		nPart = size(coords1,2);
				
		# ===================
		# Simulation results
		# ===================
		if displacement
			if fulllog
				# save variables to file
				file = matopen(string(outfilename,"_",element1,"_",nSteps,"steps_elnumber",nel,".mat"),"w")
				write(file,"HT",HT);
				write(file,"Eel0",Eel0);
				write(file,"electron_number",nel);
				write(file,"Et",Et);
				write(file, "time", time_temp)
				write(file,"coords",coord_temp)
				write(file,"coordssub",coordssub)
				write(file,"vels",vels_temp)
				write(file,"E",E_temp)
				write(file,"L",L)
				write(file,"sig1",a1)
				write(file,"eps1",eps1)
				#write(file,"sig2",a2)
				#write(file,"eps2",eps2)
				write(file,"element1",element1);
				#write(file,"element2",element2);
				write(file,"countNN",countNN_temp)
				close(file)
			end
			
			# write only geometry:
			#writeXYZ(coords1.*1E10,element1,string("afterexperiment","_eln",nel));
			writeXYZ2el(coords1.*1E10,coordssub.*1E10,element1,"C",string("afterexperiment_",element1,nPart,"C",nPartS,potential1,"_withsub","_eln",nel));
		end
		# ---------------------------------------------------------------------
		
		nel = nel + 1;
		
		println("save electron log"); 
		file = matopen(string(outfilename,"_",element1,"_",nSteps,"steps_electronlog",".mat"),"w")
		write(file,"shotel",shotelectrons_temp);
		write(file,"sputtering_list",sputteringevents);
		write(file,"displacementthreshold",displThr);
		write(file,"high_tension",HT);
		write(file,"element1",element1);
		close(file)
		
		if resetconfig
			# reset configuration:
			coords1 = zeros(nDim,nPartStart)
			coords1 = initcoords;
			# every experiment starts with same coordinates
			println("reset configuration")
		end
		
		print_with_color(:green,"files saved, ");
		toc();
	end

	print_with_color(:green,"maximum number of electrons reached\n");
