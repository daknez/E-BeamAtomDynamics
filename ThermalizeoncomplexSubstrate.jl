	
	println("-----------------------------");
	print_with_color(:green,"Thermalisation started with \n");
	println(string("Temp=",Temp,"K, nu=",nu/1E12," THz", ", Time step: ",round(dt*1E15,1)," fs"));
	println("-----------------------------");

	dt2 = dt*dt;
	time = 0; # Following simulation time
	Epot = 0;
	samplecount = 1;
	
	v_0 = [0 0 -200];	# initial (landing) velocity of the whole cluster in m/s	
	
	Ekin1 = zeros(nPart);
	Epot1 = zeros(nPart);
	countNN1 = zeros(nPart);
	
	coord_temp = zeros(nDim,nPart,round(Integer,nSteps/sampleFreq));
	vels_temp = zeros(nDim,nPart,round(Integer,nSteps/sampleFreq));
	countNN_temp = zeros(nPart,round(Integer,nSteps/sampleFreq)+1);
	time_temp = zeros(1,round(Integer,nSteps/sampleFreq));
	E_temp = zeros(2,round(Integer,nSteps/sampleFreq));
	
	# cell list parameters:
	maxnrofpairs = size(collect(combinations(1:nPart,2)),1);

	# calculate initial thermal velocities:
	velSTDT1 = sqrt(kB*Temp/m1); # Thermostat parameter - standard deviation of the velocity
	vels1 = InitVels(nPart, nDim, velSTDT1, Temp);

	vels1 = broadcast(+,vels1,v_0');	# add initial velocity of the whole cluster
	
	forces1 = zeros(size(coords1));
	
	(rhos,pairs) = calculate_rhopairs(coords1,a1, mSC1, maxd, maxnrofpairs)
	forces1 = SC_ForcePL(coords1,rhos,L,pairs,a1,eps1,mSC1,nSC1,C1, maxd);

	for stept=1:nSteps

		tempcount = 0;

		# === 1st integration step ===

		coords1 = coords1 + dt.*vels1 + 0.5.*dt2.*forces1./m1; 

        # Update velocities:
       
		vels1 = vels1 + 0.5.*dt.*forces1./m1;

		# single CPU

		# (rhos,pairs) = calculate_rhopairs(coords1,a1, mSC1, maxd, maxnrofpairs)
		# forces1 = SC_ForcePL(coords1,rhos,L,pairs,a1,eps1,mSC1,nSC1,C1, maxd);
		# (subforces,thermpartlist) = sub_LJforces(coordssub[:,1:convert(Int64,round(nPartS/3))],coords1,sigS,epsS);
		
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
		thermpartlist = thermpartlist1 | thermpartlist2 | thermpartlist3;

	
		forces1 = SC_ForcePL(coords1,rhos,L,pairs,a1,eps1,mSC1,nSC1,C1, maxd) + subforces;
		
		# Implement surface Andersen thermostat
		list = find(thermpartlist.==true)
		for part = list
			# Test for collisions with the Andersen heat bath
			if rand() < nu.*dt
				# If the particle collided with the bath, draw a new velocity out of a normal distribution
				
				# vels1[:,thermpartlist[tempcount]] = randGauss(0,velSTDT,nDim);	
				# vels1[3,thermpartlist[tempcount]] = abs(vels1[3,thermpartlist[tempcount]]);	# allow velocities only pointing away from substrate
				
				vels1[:,part] = randGauss(0,velSTDT1,nDim);	
				vels1[3,part] = abs(vels1[3,part]);			# allow only velocities away from substrate (positive z direction)
			end
			tempcount = tempcount + 1;
		end
		
		# # classical Andersen thermostat:
		# for part = 1:nPart
			# # Test for collisions with the Andersen heat bath
			# if rand() < nu.*dt
				# # If the particle collided with the bath, draw a new velocity
				# # out of a normal distribution
				# if part > nPart1
					# vels2[:,part-nPart1] = randGauss(0,velSTDT2,nDim);
				# else
					# vels1[:,part] = randGauss(0,velSTDT1,nDim);
				# end
			# end
		# end
		
        # === Second integration step ===
        # Update velocities:
        vels1 = vels1 + 0.5.*dt.*forces1./m1;

		# === Move time forward ===
        time = time + dt;
        
		# === Sample ===

		if mod(stept,sampleFreq) == 0
			
			(Ekin1, Epot1, countNN1) = calcEnergies(coords1,L,maxd,a1,eps1,m1,mSC1,nSC1,C1,vels1,rNN2);
			
			Ekin_sum = sum(Ekin1)
			Epot_sum = sum(Epot1)
		
			if mod(samplecount,100) == 0
				
				println(string("step no: ",stept))
				#println(string("ncell: ",ncell))
				println(string("Temperature of Andersen Thermostat: ",Temp," K"))
				println(string("number of atoms considered for thermostat: ",tempcount))
				println(string("elapsed time in simulation: ",time," s"))
				println(string("Ekin = ",Ekin_sum," J = ",Ekin_sum/elch," eV"))
				# Temperatur in non-periodic system 2*Ekin/3N-6 (3f translation & 3 f rotation)
				println(string("Temperature from Ekin = ",2.*Ekin_sum/(3.*nPart-6)/kB," K"))
				println(string("Epot = ",Epot_sum," J = ",Epot_sum/elch," eV"))
				#println(string("Ekin + Epot = ",sum(Ekin1+Ekin2+Epot1+Epot2)," J = ",sum(Ekin1+Ekin2+Epot1+Epot2)./elch," eV"))
				println("-------------------------------------------------------------")
			end
			# save
			time_temp[1,samplecount] = time; 
			coord_temp[1:3,:,samplecount] = coords1
			vels_temp[1:3,:,samplecount] = vels1
			E_temp[1,samplecount] = Ekin_sum;
			E_temp[2,samplecount] = Epot_sum;
			countNN_temp[:,samplecount] = countNN1
			
			samplecount = samplecount+1;
		end
	end

	# save variables to file
	file = matopen(string(outfilename,"_",element1,"_@",Temp," K_",nSteps,"steps_substrate-thermalization.mat"),"w")
	write(file, "time", time_temp)
	write(file,"coords",coord_temp)
	write(file,"vels",vels_temp)
	write(file,"E",E_temp)
	write(file,"L",L)
	write(file,"subz",subz)
	write(file,"sig1",a1)
	write(file,"eps1",eps1)
	write(file,"countNN",countNN_temp)
	write(file,"element1",element1)
	write(file,"nPartel1",nPart)
	close(file)
	
	#writeXYZ2el(coords1.*1E10,coords2.*1E10,element1,element2,string("after subtherm_steps_",nSteps));
	
	println("-----------------------------");
	print_with_color(:green,"Thermalisation finished, files written.\n");
	println("-----------------------------");