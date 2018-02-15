# ========================================================================
# This file contains various functions for molecular simulations in Julia
# by Daniel Knez, 05.02.2018
# FELMI-ZFE/University of Technology Graz
# ========================================================================

function create_rnd_cluster(nPart,density,Amass)
	
	println("creating cluster...")
	rmin = 0.25E-9; # nm
	coords = zeros(3,nPart);

	rcluster = ((3*nPart*Amass)/(4*pi*density))^(1/3);
	
	rmin2 = rmin^2;
	rmax = 1.5*rcluster;
	L = 2*rmax;
	
	accept = 1;
	count=1;

	skipped = 0;

	while count!=nPart+1
	   
		r = abs(randn()) .* 1.2 .* rcluster;
		if r>rmax
			continue;
		end
		#r= rand().* 3.*rcluster;
		theta = rand()*pi;
		phi = rand()*2*pi;
	  
		coords[:,count] = [r*sin(theta)*cos(phi);r*sin(theta)*sin(phi);r*cos(theta)];
	  
		for part=1:size(coords,2)
			
			if part == count
				continue;
			end
		  
			dr = coords[:,part] - coords[:,count];
			dr2 = dot(dr,dr);
		  
			if dr2 < rmin2
				accept = 0;
			end
		end
	  
		if accept == 1
			println(count)
			count = count + 1;
		else
			accept = 1;
			skipped = skipped + 1;
		end
		
	end
	
	return coords, L
end

function reset_atoms_over_distance(coords,rcluster,maxdistance,origin)

	nPart = size(coords,2);
	rmin2 = (1)^2;
	rmax2 = maxdistance^2;
	
	for part = 1:nPart

		distvec = coords[:,part] - origin[:,1];
		dist = dot(distvec,distvec).^0.5;
		
		if dist > maxdistance
			println(coords[:,part])
			
			accept = 0;
			while accept == 0
			
				#newcoord = [origin[1],origin[2],abs(randn())*rcluster)];
				newcoord = [origin[1]+randn()*rcluster,origin[2]+randn()*rcluster,origin[3]+abs(randn())*rcluster];
				
				if dot(newcoord-origin[:,1],newcoord-origin[:,1]) > rmax2
					continue;
				end
				for count=1:size(coords,2)
				
					if count == part
						continue;
					end
				  
					dr = coords[:,count] - newcoord;
					dr2 = dot(dr,dr);
				  
					if dr2 > rmin2
						accept = 1;
					else
						break;
					end
				end

				coords[:,part] = newcoord;
				println(coords[:,part])
				#coords = hcat(coords[:,1:part-1],coords[:,part+1:end]);		
			end
		end
	end
	return coords, nPart
end

function Morse_forces(coords,z,D,a,R,rc)
	
	# Parameters:
	nPart = size(coords,2);
	Mforces=zeros(3,nPart);
	#Epot = zeros(1,nPart);
	
	# D = 1.5; 	# in units of epsilon
	# a = 1.5;		# a=(ke/(2*D))^0.5
	# R = 1;		# position of potential minimum (covalent radius of carbon 77 pm)
	
	# Loop over all particles:
	for part = 1:nPart
		# Calculate particle-substrate distance
		dr = coords[3,part]-z;
		if dr<rc
			expf = exp(a*(R-dr));
		
			Mforces[3,part] = 2*a*D*expf*(expf-1);
			#Epot[part] = Epot[part] + (D*((1-expf)^2-1));
		end
	end
	return Mforces#, Epot
end

function sub_LJforcesAuC(coords,z,r0,D)
	
	# Parameters:
	nPart = size(coords,2);
	LJforces=zeros(3,nPart);
	#Epot = zeros(1,nPart);
	
	# D = 1.5; 	# in units of epsilon
	# a = 1.5;		# a=(ke/(2*D))^0.5
	# R = 1;		# position of potential minimum (covalent radius of carbon 77 pm)
	
	# Loop over all particles:
	for part = 1:nPart
		# Calculate particle-substrate distance
		dr = coords[3,part]-z;
		#nen2 = (dr/r0-1.2)^6;
		
		# (-1211.76+70.38 (-1.2+x/r0)^6)/(r0 (-1.2+x/r0)^13)
		LJforces[3,part] = -D*(207*(dr/r0-1.2)^6-3564)/(r0*(dr/r0-1.2)^13);
		#Epot[part] = Epot[part] + (D*((1-expf)^2-1));
	end
	return LJforces#, Epot
end

function sub_LJ93forces(coords,z,sig,eps,rhosub)
	
	# Parameters:
	nPart = size(coords,2);
	LJforces=zeros(3,nPart);
	#Epot = zeros(1,nPart);
	
	numdens = rhosub/1.99E-26;	#atomic mass carbon=12*u=1.99E-26 kg
	epsfactor = 4*pi*numdens*sig^3/3;
	
	# Loop over all particles:
	for part = 1:nPart
	
		# Calculate particle-substrate distance
		dr = coords[3,part]-z;

		# drinv2 = 1/dr^2;
		# LJforces[3,part] = epsfactor*D*rho6*(drinv2^2-drinv2^5*0.4*rho6);
		LJforces[3,part] = -3*epsfactor*eps*(5*sig^3*dr^6-2*sig^9)/(10*dr^10);
	end
	
	return LJforces#, Epot
end

function sub_LJforces(coordssub,coords,sig,eps)
	
	coordssub = coordssub./sig;
	coords = coords./sig;
	
	# Parameters:
	nPart = size(coords,2);
	nPartS = size(coordssub,2);
	LJforces = zeros(3,nPart);
	Epot=0;
	maxd = 9;#16;
	thermpartlist = falses(nPart);
	rNN2 = 1.44;#1.3^2;
	
	# Loop over all particles:
	for part = 1:nPart
		for partS = 1:nPartS
		
			dr = coordssub[:,partS] - coords[:,part];
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]
		
			if dr2<maxd	

				invDr2 = 1/dr2; # sigma/r^2
				invDr6 = invDr2^3;
				forceFact = invDr2^4 * (invDr6 - 0.5);

				LJforces[:,part] = LJforces[:,part] - dr*forceFact;
				
				if dr2<rNN2;
					thermpartlist[part] = true;
				end
				#Epot = Epot + invDr6 * (invDr6-1);
			end
		end
	end
	LJforces = 48.*LJforces.*(eps./sig);	# convert forces to SI - units
	#LJforces = 48.*eps.*LJforces;
	# delete double elements from thermalisation list:
	#thermpartlist = [unique(thermpartlist);zeros(nPart-length(thermpartlist))];
	
	return LJforces,thermpartlist#, Epot
end

function sub_potentialforces(coords1,Lx,Ly,Lz,forcex,forcey,forcez,subz)
	nDim = size(coords1,1);
	nPart = size(coords1,2);
	forces = zeros(nDim,nPart);
	# sampling distance in each direction:
	psx = Lx/(size(forcex,1)-1);
	psy = Ly/(size(forcey,1)-1);
	psz = Lz/(size(forcez,1)-1);
	x = collect(linspace(0,Lx,size(forcex,1)));
	y = collect(linspace(0,Ly,size(forcey,1)));
	z = collect(linspace(0,Lz,size(forcez,1)));
	sig = 0.3e-9;
	
	# Loop over all particles:
	for part = 1:nPart
		# Calculate particle-substrate distance
		#dr = [coords1[1,part],coords1[2,part],coords1[3,part]-subz];
		dr = coords1[:,part];
		
		indexx = find(dr[1].<x.<=dr[1]+psx);
		indexy = find(dr[2].<y.<=dr[2]+psy);
		indexz = find(dr[3].<z.<=dr[3]+psz);
		forces[1,part] = -forcex[indexx,indexy,indexz][1];
		forces[2,part] = -forcey[indexx,indexy,indexz][1];
		#forces[3,part] = -48*potsurf[indexx,indexy][1]*((sig/dr[3])^12-0.5*(sig/dr[3])^6);
		#forces[3,part] = -24*potsurf[indexx,indexy][1]/sig*(2*(sig/dr[3])^13-(sig/dr[3])^7);
		forces[3,part] = forcez[indexx,indexy,indexz][1];
	end

return forces
end

function initCubicGrid(nPart,density)

	# Initialize with zeroes
	coords = zeros(3,nPart);

	# Get the cooresponding box size
	L = (nPart/density)^(1.0/3);

	# Find the lowest perfect cube greater than or equal to the number of
	# particles
	nCube = 2;
	
	while (nCube^3 < nPart)
		nCube = nCube + 1;
	end
	
	# Start positioning - use a 3D index for counting the spots
	index = [0,0,0]';
	
	# Assign particle positions
	for part=1:nPart
		# Set coordinate
		coords[:,part] = (index+[0.5,0.5,0.5]')*(L/nCube);
		
		# Advance the index
		index[1] = index[1] + 1;
		if (index[1] == nCube) 
			index[1] = 0;
			index[2] = index[2] + 1;
			if (index[2] == nCube)
				index[2] = 0;
				index[3] = index[3] + 1;
			end
		end
	end
	
	coords = coords + L;
	L=3*L;
return coords,L
end

function cell_list(pos,ncel,rcel) #generating cell_list 

	# pos... coordinates of each particle (nDim,nPart)
	# ncel... number of cells per dimension 
	# rcel... size of a cell 
	# n... label of particles 
	# cellx,celly,cellz... id of the cell for certain atom 
	# link(nPart) 
	# cell(1:ncel,1:ncel,1:ncel) 
	 
	nDim = size(pos,1);
	nPart = size(pos,2);
	
	# initialize link- and cell-list
	link = round(Integer,zeros(nPart));

	celll = zeros(ncel,ncel,ncel);
	 
	# avoid negative indices:
	pos = broadcast(+,pos,abs(minimum(pos)));
	
	#set up the cell list 
	for n = 1:nPart #scan all atoms 
		#locate the atom to certain cell 
		cellx = round(Integer,pos[1,n]./rcel)+1; 
		celly = round(Integer,pos[2,n]./rcel)+1; 
		cellz = round(Integer,pos[3,n]./rcel)+1; 
		#insert the atom into the link list of the corresponding cell 
		
		link[n] = celll[cellx,celly,cellz];	# link list: in which cell can each atom be found
		celll[cellx,celly,cellz] = n; 		# cell list: contains the head of the corresponding link list
	end
	 
	return link,celll
end


function LJForcesPar(coords,L)
	nPart = size(coords,2);
	forces = zeros(size(coords));
	parad1 = @spawn LJ_Force(coords[:,1:round(Integer,nPart./2)],L);
	parad2 = @spawn LJ_Force(coords[:,round(Integer,nPart./2)+1:end],L);
	parad3 = @spawn LJ_ForceAB(coords,L);
		
	forces[:,1:round(Integer,nPart./2)] = fetch(parad1);
	forces[:,round(Integer,nPart./2)+1:end] = fetch(parad2);
	forces = forces + fetch(parad3);
	
return forces
end

function InitVels(nPart, nDim, velSTDT, TempT)

	
	# Set initial velocities with random numbers
	vels1 = zeros(nDim,nPart);
	for part=1:nPart
		vels1[:,part] = randGauss(0,velSTDT,nDim);
	end

	# Set initial momentum to zero
	totV = sum(vels1,2)/nPart; # Center of mass velocity
	for dim=1:nDim
		vels1[dim,:] = vels1[dim,:] - totV[dim];    # Fix any center-of-mass drift
	end

	# # Set initial kinetic energy to nDim*KbT/2
	# totV2 = (vels1[1,:]*vels1[1,:]'+vels1[2,:]*vels1[2,:]'+vels1[3,:]*vels1[3,:]')/nPart; # mean sqared velocity
	# velScale = sqrt(3*nPart.*TempT./totV2);   # Velocity scaling factor
	# vels1 = vels1.*velScale;    


return vels1
end

function calculate_rho(coords, m, a, maxd)

	nPart = size(coords,2);
	rho = zeros(nPart);
	
	for partA = 1:nPart-1
		for partB = (partA+1):nPart
			
			dr = coords[:,partA] - coords[:,partB];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]
			
			if dr2<maxd
				sumrho = (a/(dr2^0.5))^m;
				rho[partA] = rho[partA] + sumrho;
				rho[partB] = rho[partB] + sumrho;
			end
		end
	end
	
	rhos = rho.^(-0.5);
	
	return rhos
end

# calculate rho for Sutton-Chen potential and return a pair list:
function calculate_rhopairs(x,a, m, maxd, maxnrpair)

	nPart = size(x,2);
	rho = zeros(nPart);
	pairs = round(Integer,zeros(maxnrpair, 2));
	pairlen = 1;
	#ps = LUTrho[2,1]-LUTrho[1,1];
	
	for partA = 1:nPart-1
		for partB = (partA+1):nPart
			
			dr = x[:,partA] - x[:,partB];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]
			
			if dr2<maxd
				
				#index = find(LUTrho[:,1].<dr2.<=LUTrho[:,1]+ps);
				#sumrho = LUTrho[index,2][1];
				sumrho = (a/(dr2^0.5))^m;
				
				rho[partA] = rho[partA] + sumrho;
				rho[partB] = rho[partB] + sumrho;
				pairs[pairlen,:] = [partA,partB];
				pairlen = pairlen + 1;
			end
		end
	end
	
	rhos = rho.^(-0.5);
	pairs = pairs[1:pairlen-1,:];

	return rhos, pairs
	#return rho, pairs
end

# Calculate rho for Sutton-Chen potential with 2 elements
function calculate_rhopairs2el(coords1,a1,m1,coords2,a2,m2,a12,m12,maxd,maxnrpair)

	x = [coords1 coords2];

	nPart1 = size(coords1,2);
	nPart2 = size(coords2,2);
	nPart = size([coords1 coords2],2);
	
	rho1 = zeros(nPart1);
	rho2 = zeros(nPart2);
	rho = zeros(nPart);
	
	pairs = round(Integer,zeros(maxnrpair,2));

	pairlen = 1;
	
	for partA = 1:nPart-1
		for partB = (partA+1):nPart
			
			dr = x[:,partA] - x[:,partB];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]
			
			if dr2<maxd
				if partA < nPart1
					if partB < nPart1
						sumrho = (a1/(dr2^0.5))^m1;
						rho[partA] = rho[partA] + sumrho;
						rho[partB] = rho[partB] + sumrho;
					else
						sumrho = (a12/(dr2^0.5))^m12;
						rho[partA] = rho[partA] + sumrho;
						rho[partB] = rho[partB] + sumrho;
					end
				else
					if partB < nPart1
						sumrho = (a12/(dr2^0.5))^m12;
						rho[partA] = rho[partA] + sumrho;
						rho[partB] = rho[partB] + sumrho;
					else
						sumrho = (a2/(dr2^0.5))^m2;
						rho[partA] = rho[partA] + sumrho;
						rho[partB] = rho[partB] + sumrho;
					end
				end
				pairs[pairlen,:] = [partA,partB];
				pairlen = pairlen + 1;
			end
		end
	end
	
	# rhos1 = rho1.^(-0.5);
	# rhos2 = rho2.^(-0.5);
	rhos = rho.^(-0.5);
	pairs = pairs[1:pairlen-1,:];

	return rhos,pairs
end

function SC_ForcePL2el(coords1,coords2,rhos,L,pairs,a1,a2,a12,eps1,eps2,eps12,m1,m2,m12,n1,n2,n12,C1,C2,maxd)	

	nPart1 = size(coords1,2);

	X = [coords1 coords2];
	nPart = size(X,2);
	forces = zeros(size(X));
	
	# Loop over all particle pairs
	pairslen = size(pairs,1);
	
	for pair = 1:pairslen
			# Calculate particle-particle distance
			dr = X[:,pairs[pair,1]] - X[:,pairs[pair,2]];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]
			
			if pairs[pair,1] < nPart1
				if pairs[pair,2] < nPart1
					drs = a1/(dr2^0.5);
					forceFact = (eps1*n1*drs^n1-0.5*eps1*C1*m1*drs^m1*(rhos[pairs[pair,1]]+rhos[pairs[pair,2]]))/dr2;
				else
					drs = a12/(dr2^0.5);
					forceFact = (eps12*n12*drs^n12-0.5*eps1*C1*m12*drs^m12*(rhos[pairs[pair,1]]+rhos[pairs[pair,2]]))/dr2;
				end
			else
				if pairs[pair,2] < nPart1
					drs = a12/(dr2^0.5);
					forceFact = (eps12*n12*drs^n12-0.5*eps2*C2*m12*drs^m12*(rhos[pairs[pair,1]]+rhos[pairs[pair,2]]))/dr2;
				else
					drs = a2/(dr2^0.5);
					forceFact = (eps2*n2*drs^n2-0.5*eps2*C2*m2*drs^m2*(rhos[pairs[pair,1]]+rhos[pairs[pair,2]]))/dr2;
				end
			end
			
			forces[:,pairs[pair,1]] = forces[:,pairs[pair,1]] + dr*forceFact;
			forces[:,pairs[pair,2]] = forces[:,pairs[pair,2]] - dr*forceFact;
	end
	
	return forces
end

function SC_ForcePL(coords,rhosqrtr,L,pairs,a,epsilon,m,n,C, maxd)	

	forces = zeros(size(coords));
	nPart = size(coords,2);
	
	# calculate all possible pairs:
	#comb = collect(combinations(1:nPart,2));
	#comb = combinations(1:nPart,2);
	# rhosqrtr = zeros(nPart);
	# rhosqrtr = calculate_rho(coords, m, a, maxd);
	
	# Loop over all particle pairs
	pairslen = size(pairs,1);
	
	for pair = 1:pairslen
			# Calculate particle-particle distance
			dr = coords[:,pairs[pair,1]] - coords[:,pairs[pair,2]];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]

			drs = a/dr2^0.5;
				
			forceFact = (n*drs^n-0.5*C*m*drs^m*(rhosqrtr[pairs[pair,1]]+rhosqrtr[pairs[pair,2]]))/dr2;
			
			forces[:,pairs[pair,1]] = forces[:,pairs[pair,1]] + dr*forceFact;
			forces[:,pairs[pair,2]] = forces[:,pairs[pair,2]] - dr*forceFact;
	end
	
	forces = epsilon .* forces;
	return forces
end

function SC_Force(coords,L,a,epsilon,m,n,C, maxd)	

	forces = zeros(size(coords));
	nPart = size(coords,2);
	
	# calculate all possible pairs:
	#comb = collect(combinations(1:nPart,2));
	#comb = combinations(1:nPart,2);
	rhosqrtr = zeros(nPart);
	rhosqrtr = calculate_rho(coords, m, a, maxd);
	forcerep = 0;
	forceattr = 0;
	# Loop over all particle pairs
	for partA = 1:nPart-1
		for partB = (partA+1):nPart
			# Calculate particle-particle distance
			dr = coords[:,partA] - coords[:,partB];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]

			if dr2<maxd	
				
				drs = a/dr2^0.5;
				
				forceFact = (n*drs^n-0.5*C*m*drs^m*(rhosqrtr[partA]+rhosqrtr[partB]))/dr2;
				#forcerep = (n*drs^n)/dr2;
				#forceattr = -0.5*C*m*drs^m*(rhosqrtr[partA]+rhosqrtr[partB])/dr2;
				forces[:,partA] = forces[:,partA] + dr*forceFact;
				forces[:,partB] = forces[:,partB] - dr*forceFact;

			end

		end
	end
	#println(string(epsilon*forcerep,", ",epsilon*forceattr));
	forces = epsilon .* forces;
	return forces
end

function calcEnergies(coords,L,maxd,a,epsilon,mass,m,n,C,vels,NNthres)
	
	nPart = size(coords,2);
	#NNthres = 0.723; # 0.85.^2
	
	countNN = zeros(nPart);
	Ekin = zeros(nPart);
	Epot = zeros(nPart);

	rhosqrtr = zeros(nPart);
	rhosqrtr = calculate_rho(coords, m, a, maxd);
	
	# Loop over all particle pairs
	for partA = 1:nPart-1
		for partB = (partA+1):nPart
			# Calculate particle-particle distance
			dr = coords[:,partA] - coords[:,partB];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]

			if dr2<maxd	
				
				drs = a/(dr2^0.5);
				
				Epot[partA] = Epot[partA] + drs^n - C* rhosqrtr[partA];
				Epot[partB] = Epot[partB] + drs^n - C* rhosqrtr[partB];
				
				if dr2<NNthres
					countNN[partA] = countNN[partA] + 1;
					countNN[partB] = countNN[partB] + 1;
				end
			end

		end
	end
			
	for i = 1:nPart
		Ekin[i] = dot(vels[:,i],vels[:,i])*0.5*mass;
	end
	Epot = Epot * epsilon;
	
	return Ekin, Epot, countNN
end



function calcEpotSC(coords,L,maxd,m,n,C)
	
	a = 1;
	epsilon = 1;

	Epot = 0;
	#Epot2 = 0;
	nPart = size(coords,2);
	
	rho = zeros(nPart);
	rho = calculate_rho(coords, m, a, maxd);
	
	# Loop over all particle pairs
	
	for partA = 1:nPart
		for partB = 1:nPart
			if partA == partB
				continue
			end
			
			# Calculate particle-particle distance
			dr = coords[:,partA] - coords[:,partB];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]

			if dr2<maxd	
				
				drs = 1/(dr2^0.5);
				
				Epot = Epot + 0.5*(a*drs)^n;
				#Epot2 = Epot2 + (a*drs^n - C* rhosqrtr[partA]).^2;
			end
		end
		Epot = Epot - C*rho[partA].^0.5;
	end
	
	return Epot#, Epot2
end

function calcdrspairs(coords,pairs)

	pairslen = size(pairs,1);
	drspairs = zeros(pairslen);
	
	for pair = 1:pairslen
	
			# Calculate particle-particle distance
			dr = coords[:,pairs[pair,1]] - coords[:,pairs[pair,2]];
			
			dr2 = dot(dr,dr);
			drs = 1/(dr2^0.5);	
			drspairs[pair] = drs;
	end
	
	return drspairs
end

function create_neighbourlist(x, rc, srgdat, frc, maxnebr)

	x = x';

	nx = size(x,1);
	maxd = rc^2;
	
	region = srgdat[1,4:6];
	(celllist, ncell, cell_offset, noffset, id_cell) = make_celllist_cube(region, 1.2*rc, broadcast(-,x,srgdat[1,1:3]), nx, frc);

	# init.
	# maxnebr = 200 * nx;
	maxnebr = 2*maxnebr;
	nebrTab = zeros(maxnebr, 1);
	nebrLen = 0;

	# creating neighborlist
	for k = 1:ncell[3]
		for j = 1:ncell[2]
			for i = 1:ncell[1]
				# 1st cell
				m1 = round(Integer,(k-1)*ncell[2]*ncell[1] + (j-1)*ncell[1] + i);
				for l = 1:noffset
					# 2nd cell
					m2v = [i j k] + cell_offset[l, :];
					p1  = celllist[ m1 + nx];
					if all( ncell - m2v .>= 0) && all( m2v .> 0)
						m2  = round(Integer,(m2v[3]-1)*ncell[2]*ncell[1] + (m2v[2]-1)*ncell[1] + m2v[1]);
						p2  = celllist[m2 + nx];
						while p1 != -1
							while p2 != -1
								if m1 != m2 || p1 > p2
								
									dr = vec(x[p1, :] - x[p2, :]);
									dr2 = dot(dr,dr);
									if dr2 <= maxd
									    nebrLen = nebrLen + 2;
										if nebrLen > maxnebr
											println(string("too many neighbours (",nebrLen," > ",maxnebr,")"));
											println();
											return;
										end
										nebrTab[nebrLen - 1] = p1;
										nebrTab[nebrLen]     = p2;
									end
								end
								p2 = celllist[p2];
							end
							p2 = celllist[m2 + nx];
							p1 = celllist[p1];
						end
					end
				end
			end
		end
	end
	
	nebrTab = round(Integer,nebrTab[1:nebrLen]);

	nebrTab = [nebrTab[1:2:end] nebrTab[2:2:end]];
	nebrLen = nebrLen / 2;

	return nebrTab, nebrLen, ncell, id_cell
end

function SetMatParameters(element1,potential)

if potential == "SC"
	if element1 == "Ni"
		# eps1 = epsNi;
		# sig1 = sigNi;
		m = 58.6934*u;
		
		# Sutton-Chen (aus Todd.1993):
		mSC=6;
		nSC=9;
		C=39.432;
		eps = 182.29 * kB;
		sig = 0.352E-9;
		
	elseif element1 == "Au"
		m = 3.27E-25;
		
		# Sutton-Chen (aus Sutton & Chen 1990):
		mSC=8;
		nSC=10;
		C=34.408;
		eps = 0.012793 * elch;
		sig = 0.408E-9;
	elseif element1 == "Ag"
		m = 1.791E-25;
		
		# Sutton-Chen (aus Sutton & Chen 1990):
		mSC=6;
		nSC=12;
		C=144.41;
		eps = 0.0025415 * elch;
		sig = 0.409E-9;
	elseif element1 == "Cu"
		m = 63.546*u;;
		
		# Sutton-Chen (aus Sutton & Chen 1990):
		mSC=6;
		nSC=12;
		C=39.432;
		eps = 1.2382E-2 * elch;
		sig = 0.361E-9;
	elseif element1 == "Pt"
		m = 195.08*u;
		
		# quantum Sutton-Chen 
		mSC = 8;
		nSC = 10;
		C = 71.336;
		eps = 230.14*kB;
		sig = 0.392E-9;
	end
elseif potential == "LJ"

	# Lennard-Jones potentials:
	if element1 == "Au"
		eps =  0.229*elch;		# in J nach Heinz et al. "Accurate Simulation of Surfaces and Interfaces of Face-Centered Cubic Metals Using 12−6 and 9−6 Lennard-Jones Potentials"
		sig = 0.2951E-9; 		# sigma in m
		m = 3.27E-25;     		# gold mass in kg
	
	elseif element1 == "Ag"
		#epsAg = 0.345.*elch;		# nach Guan et al. "MD simulations of Ag ﬁlm growth using the Lennard-Jones potential
		#sigAg = 0.2644E-9;		# Guan et al.
		eps = 0.1977*elch;		# nach Heinz et al.
		sig = 0.2955E-9;		# Heinz et al
		m = 107.86*u;
	
	elseif element1 == "Ni"
		eps = 0.245*elch;		# Heinz et al.
		sig = 0.2552E-9;		# Heinz et al.
		m = 58.6934*u;
	end
	mSC = nSC = C = 0;
	
elseif potential == "QSC"	
#T. Çagin, Y. Qi, H. Li, Y. Kimura, H. Ikeda, W.L. Johnson, W.A. Goddard III, The quantum Sutton-Chen many-body potential for properties of fcc metals, MRS Symposium Ser. 554 (1999) 43. 
	if element1 == "Ni"
		m = 58.6934*u;
		
		# quantum Sutton-Chen
		mSC = 5;
		nSC = 10;
		C = 84.745;
		eps = 7.3767E-3 * elch;
		sig = 0.35157E-9;
		
	elseif element1 == "Au"
		m = 3.27E-25;
		
		#quantum Sutton-Chen:
		mSC = 8;
		nSC = 11;
		C = 53.581;
		eps = 7.8052E-3 * elch;
		sig = 0.40651E-9;
	elseif element1 == "Ag"
		m = 1.791E-25;
		
		mSC = 6;
		nSC = 11;
		C = 96.524;
		eps = 3.9450E-3 * elch;
		sig = 0.40691E-9;
	
	elseif element1 == "Cu"
		m = 63.546*u;
		
		# quantum Sutton-Chen
		mSC = 5;
		nSC = 10;
		C = 84.843;
		eps = 5.7921E-3 * elch;
		sig = 0.36030E-9;
	elseif element1 == "Pt"
		m = 195.08*u;
		
		# quantum Sutton-Chen 
		mSC = 7;
		nSC = 11;
		C = 71.336;
		eps = 9.7894E-3 * elch;
		sig = 0.39163E-9;
	end
end
	
return sig, eps, mSC, nSC, C, m
end

#----------------------------------------------------------------
# write Output files
#----------------------------------------------------------------

function writeXYZ(X,element,filename)
	
	nPart = size(X,2);
	#X =  coords1.*sig1.*1E10;
	file = open(string(filename,".xyz"),"w")
	write(file,string(nPart,"\n\n"));
	for ii = 1:nPart;
		#write(file,X[1,ii],X[2,ii],X[3,ii],"Au","\n");
		write(file,string(element," ",X[1,ii]," ",X[2,ii]," ",X[3,ii],"\n"));
	end
	close(file);
	
	println(string(filename,".xyz"," written!"));
end

function writeXYZ2el(X1,X2,element1,element2,filename)
	
	nPart1 = size(X1,2);
	nPart2 = size(X2,2);

	file = open(string(filename,".xyz"),"w")
	write(file,string(nPart1+nPart2,"\n\n"));
	for ii = 1:nPart1;
		write(file,string(element1," ",X1[1,ii]," ",X1[2,ii]," ",X1[3,ii],"\n"));
	end
	for ii = 1:nPart2;
		write(file,string(element2," ",X2[1,ii]," ",X2[2,ii]," ",X2[3,ii],"\n"));
	end
	close(file);
	
	println(string(filename,".xyz"," written!"));
end

#---------------------------------------------------------------------------------
# ---------------------------- MC related Functions ------------------------------
#---------------------------------------------------------------------------------

function calcdeltV_SC(coords,trialPos,mSC,nSC,C,maxd,part)
		
		a = 1;
	
		deltaE = 0;
		deltaErep = 0;
		repENew = 0;
		repEOld = 0;
    
        # Get the number of particles
        nPart = size(coords,2);
		
		coords_New = zeros(size(coords));
		coords_New = deepcopy(coords);
		coords_New[:,part] = trialPos;
        
		
		# rho_Old = zeros(nPart);
		# rho_New = zeros(nPart);

		rho_Old = 0;
		rho_New = 0;
		
		# rho_Old = calculate_rho(coords, mSC, a, maxd);
		# rho_New = calculate_rho(coords_New, mSC, a, maxd);
		rho_Old = calculate_rhoMC(coords, mSC, a, maxd,part);
		rho_New = calculate_rhoMC(coords_New, mSC, a, maxd,part);
		
		# println(C*rho_Old^0.5)
		# println(C*rho_New^0.5)
		
        for otherPart = 1:nPart
            
            # Make sure to skip particle 'part' so that we don't calculate self
            # interaction
            if otherPart == part
                continue
            end
                
            # Calculate particle-particle distance for both the old and new
            # configurations
            drNew = coords[:,otherPart] - coords_New[:,part];
            drOld = coords[:,otherPart] - coords[:,part];
            
            # Get the distance squared
            dr2_New = dot(drNew,drNew);
            dr2_Old = dot(drOld,drOld);
			
            # ENew = (a*drs_New)^nSC - C* rhosqrtr_New[part];
			# EOld = (a*drs_Old)^nSC - C* rhosqrtr_Old[part];
			# ENew = (a*drs_New)^nSC - C* rhosqrtr_New;
			# EOld = (a*drs_Old)^nSC - C* rhosqrtr_Old;
			# ENew = (a/(dr2_New^0.5))^nSC - C*(rho_New)^0.5;
			# EOld = (a/(dr2_Old^0.5))^nSC - C*(rho_Old)^0.5;
			
			#ENew = 0.5*(a/(dr2_New^0.5))^nSC;
			#EOld = 0.5*(a/(dr2_Old^0.5))^nSC;
			
			# repENew = repENew + 0.5*(a/sqrt(dr2_New))^nSC;
			# repEOld = repEOld + 0.5*(a/sqrt(dr2_Old))^nSC;
			
			repENew = repENew + 0.5*(a/sqrt(dr2_New))^nSC;
			repEOld = repEOld + 0.5*(a/sqrt(dr2_Old))^nSC;
			
        end
		ENew = repENew - C*(rho_New)^0.5
		EOld = repEOld - C*(rho_Old)^0.5
		
		deltaE = ENew - EOld;
		
		return deltaE
end

function calculate_rhoMC(coords, mSC, a, maxd, part)

	nPart = size(coords,2);
	rho = 0;
	
	for partB = 1:nPart
	
		if partB == part
			continue
		end
			
		dr = coords[:,partB] - coords[:,part];
	
		dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]
		
		if dr2<maxd
			sumrho = (a/(dr2^0.5))^mSC;
			rho = rho + sumrho;
		end
	end
	
	return rho
end

function sub_LJdeltaEAuC(oldPos,trialPos,z,r0,D,part)

	deltaE = 0;
	
	drOld = oldPos[3] - z;
	drNew = trialPos[3] - z;
	
	# Get the distance squared
	dr_New = dot(drNew,drNew)^0.5;
	dr_Old = dot(drOld,drOld)^0.5;
	
	ENew = D * (297/(dr_New/r0-1.2)^12 - 34.5/(dr_New/r0-1.2)^6);
	EOld = D * (297/(dr_Old/r0-1.2)^12 - 34.5/(dr_Old/r0-1.2)^6);
	
	deltaE = ENew - EOld;
	
	return deltaE
end

 function LJ_EnergyChange(coords,trialPos,part,L)
    
        deltaE = 0;
    
        # Get the number of particles
        nPart = size(coords,2);
        
        # Loop over all particles and calculate interaction with particle
        # 'part'.
        for otherPart = 1:nPart
            
            # Make sure to skip particle 'part' so that we don't calculate self
            # interaction
            if otherPart == part
                continue
            end
                
            # Calculate particle-particle distance for both the old and new
            # configurations
            drNew = coords[:,otherPart] - trialPos;
            drOld = coords[:,otherPart] - coords[:,part];
            
            # Get the distance squared
            dr2_New = dot(drNew,drNew);
            dr2_Old = dot(drOld,drOld);
    
            # Lennard-Jones potential:
            # U(r) = 4*epsilon* [(sigma/r)^12 - (sigma/r)^6]
            #
            # Here, we set sigma = 1, epsilon = 1 (reduced distance and
            # energy units). Therefore:
            #
            # U(r) = 4 * [(1/r)^12 - (1/r)^6]
            # 
            # For efficiency, we will multiply by 4 only after summing
            # up all the energies.
                    
            invDr6_New = 1.0/(dr2_New^3); # 1/r^6
            invDr6_Old = 1.0/(dr2_Old^3); # 1/r^6
            
            # Calculate the potential energy
            eNew = (invDr6_New * (invDr6_New - 1));
            eOld = (invDr6_Old * (invDr6_Old - 1));
            
            deltaE = deltaE + eNew - eOld;
        end
        
        # Multiply energy by 4
        deltaE = deltaE*4;
		
		return deltaE
    end

function calcNN(coords,NNthres)
	
	nPart = size(coords,2);
	
	countNN = zeros(nPart);

	# Loop over all particle pairs
	for partA = 1:nPart-1
		for partB = (partA+1):nPart
			# Calculate particle-particle distance
			dr = coords[:,partA] - coords[:,partB];
			# Fix according to periodic boundary conditions
			# dr = distPBC3D(dr,L);
			
			dr2 = dot(dr,dr) # = dr[1].*dr[1]+dr[2].*dr[2]+dr[3].*dr[3]

			if dr2<NNthres
				countNN[partA] = countNN[partA] + 1;
				countNN[partB] = countNN[partB] + 1;
			end

		end
	end
			
	return countNN
end

function randGauss(mu,sigma,nDim)
	# Generate normally distributed random numbers
	randNums = randn(nDim,1);
	
	# Shift to match given mean and std
	randNums = mu + randNums * sigma;
	return randNums
end

function rndPowerLaw(min,max,a,b)
	
	numberfound = false;
	randnum = 0;
	while numberfound != true
		randnum = (-rand()*(b+1)/a).^(1/(b+1));
		if randnum > min && randnum < max 
			numberfound = true;
		end
	end
	
	return randnum
end
function MC_Cluster_relaxation(coords1,L,element1,potential,trialsperatom,sampleFreq,startTemp,outfilename)
	
	(sig1, eps1, mSC, nSC, C, m) = SetMatParameters(element1,potential);
	
	nPart = size(coords1,2);
	nDim = size(coords1,1);
	coords1=coords1./sig1;
	#outfilename = string("MCrelaxation_",nPart,element1);
	
	t1=t2=t3=0;

	mass = 1.0;
	sqrt2 = sqrt(2)/2;
	rNN = sqrt2*1.2;
	rNN2 = rNN^2;
	rc = 2.5 * rNN; 	# distances are in units of sigma which is a in fcc in SC-units (rc = 2.5*rNN, sqrt(2)/2 = APF for fcc)
	#maxd = rc^2;
	maxd = 100;
	
	#trialsperatom = 3000;
	Ntrials = trialsperatom*nPart;
	#sampleFreq = 1000;

	Temp = startTemp / (eps1/kB);
	beta = 1.0/Temp;
	# parameters for simulated annealing
	tempstep = 0.2 / (eps1/kB);
	tempchangeinterval = nPart;	
	templowlimit = 0.4 / (eps1/kB);
	drmax = 0.1*rNN
	
	Naccepted = 0;
	Nskipped = 0;
	samplecount = 1;
	acceptratio = 0;
	atom = 1;

	countNN = calcNN(coords1,rNN2);
	
	coord_temp = zeros(size(coords1,1),size(coords1,2),round(Integer,Ntrials/sampleFreq));
	energy_temp = zeros(round(Integer,Ntrials/sampleFreq));
	countNN_temp = zeros(nPart,round(Integer,Ntrials/sampleFreq));
	steps=1;

	writeXYZ(coords1.*sig1.*1E10,element1,string(outfilename,"_preMC"));
	println(Ntrials/nPart);
	energy = calcEpotSC(coords1,L,maxd,mSC,nSC,C);
	tic();

	for steps=1:Ntrials
		
		# simulated annealing:
		if mod(steps,tempchangeinterval) == 0 && Temp > templowlimit+tempstep
			Temp = Temp - tempstep;
			beta = 1/Temp;
		elseif Temp-tempstep < templowlimit
			Temp = templowlimit;
			beta = 1/Temp;
		end
			
			
		#for atom=1:nPart
		# only consider surface atoms (surface diffusion):
			atom = round(Integer,floor(rand()*nPart)+1);
			
			# if countNN[atom] > 5
				# Nskipped = Nskipped + 1;
				# continue;
			# end
			
			#rnew = (2.*rand(nDim)-[1;1;1]).*drmax;
			rnew = drmax.*(rand(nDim)-0.5);
			trialpos = coords1[:,atom] + rnew;
			
			coords1_trial = deepcopy(coords1);
			coords1_trial[:,atom] = trialpos;
			# Vold = calcV_SC(coords1,mSC,nSC,C,maxd,atom);
			# Vnew = calcV_SC(coords1_trial,mSC,nSC,C,maxd,atom);
			# deltV = Vnew - Vold;
			
			#if potential == "SC"
			deltV = calcdeltV_SC(coords1,trialpos,mSC,nSC,C,maxd,atom);
			#elseif potential == "LJ"
			#	deltV = LJ_EnergyChange(coords1,trialpos,atom,L)
			#end
			
			# LJ Au-aC Interaction:
			#deltV = deltV + sub_LJdeltaEAuC(coords1[:,atom],trialpos,subz,r0,DM,atom);
			
			deltVb = beta * (deltV);

			if deltVb < 50
				if deltVb < 0.0
					coords1[:,atom] = trialpos;
					energy = energy + deltV;
					Naccepted = Naccepted + 1;
				elseif exp(-deltVb) > rand()
					coords1[:,atom] = trialpos;
					energy = energy + deltV;
					Naccepted = Naccepted + 1;
				end
			end
		#end
		if mod(steps,sampleFreq) == 0
			#acceptratio = Naccepted/(steps-Nskipped);
			acceptratio = Naccepted/sampleFreq;
			if acceptratio > 0.6
				drmax = drmax * 1.05;
				println(string("maximum displacement changed to ",drmax));
			elseif acceptratio < 0.5
				drmax = drmax * 0.95;
				println(string("maximum displacement changed to ",drmax));
			end
			Naccepted = Nskipped = 0;
			
			coord_temp[1:3,:,samplecount] = coords1.*sig1;
			energy_temp[samplecount] = energy.*eps1;
			#countNN = calcNN(coords1,rNN2)
			samplecount = samplecount + 1;
			println(string("Step: ",steps," of ",Ntrials," Temperature: ",Temp*(eps1/kB)," ",acceptratio*100," %"))
			t1 = toc();
			tic();
			println(string((Ntrials-steps)*t1/3600/sampleFreq," hours remaining..."))
			# if mod(steps,sampleFreq*100) == 0
				# writeXYZ(coords1.*sig1.*1E10,element1,string(outfilename,"_",element1,"_",Ntrials,"_",steps));
			# end
			
		end
		
	end

	writeXYZ(coords1.*sig1.*1E10,element1,string(outfilename,"_MCresult_trialsteps",Ntrials));

	println(energy_temp[1])
	println(energy_temp[end])
	#println(string(steps," ",atom," ",deltV," ",acceptratio*100," %"))

	# save variables to file
	file = matopen(string(outfilename,"_",element1,"_",Ntrials,".mat"),"w")
	write(file,"coords1",coord_temp)
	write(file,"energy",energy_temp)
	write(file,"countNN",countNN_temp)
	write(file,"L",L)
	#write(file,"subz",subz.*sig1)
	write(file,"sig1",sig1)
	write(file,"eps1",eps1)
	close(file)
	
return coords1.*sig1
end

function create_aCsubstrate(L,d,density)

	rmin = 0.145E-9; # nm
	# L=4;
	# d=0.5;
	range = [L ; L ; d]; # m
	# density = 2200; # kg/m³
	Amass = 12.011*u; # kg

	rmin2 = rmin^2;
	nPart = round(Integer,range[1]*range[2]*range[3]*density./Amass);
	accept = 1;
	count=1;
	coords = zeros(3,nPart);

	skipped = 0;

	while count!=nPart+1
	   
	  coords[:,count] = rand(3,1).*range;
	  
	  # thickness normal distributed:
	  # coords(1:2,count) = rand(2,1).*range(1:2);
	  # coords(3,count) = randn().*range(3)*0.5;
	  
	  for part=1:size(coords,2)
			
		  if count == part
			  continue;
		  end
		  
		  dr = coords[:,count] - coords[:,part];
		  dr2 = dot(dr,dr);
		  
		  if dr2 < rmin2
			  accept = 0;
		  end
	  end
	  
	  if accept == 1
		count = count + 1
	  else
		  accept = 1;
		  skipped = skipped + 1;
	  end
		
	end

return coords
end

function calcscatteringvelvector(m,E,theta)

	Et = 2*Eel*(Eel+2*Eel0)/(c^2*m) * sin(theta/2)^2; 		# transferred Energy (electron to atom) in J

	ve = c*(1-Eel0^2/(Eel0+Eel-Et)^2)^0.5;	# speed of electron after scattering in LJ-units
	ve0 = c*(1-Eel0^2/(Eel0+Eel)^2)^0.5;	# speed of electron before scattering

	vatom = (2*Et/m)^0.5; 				# speed of hit atom after scattering 
	
	psi = asin(mel*ve/(1-(ve/c)^2)^0.5/(m*vatom)*sin(theta));	# from momentum conservation (p_el*sin(theta) = p_atom *sin(psi) & relativistic p_el and classic pa=m*vatom
	phi = rand()*2*pi;

	v_sc = [cos(phi)*sin(psi); sin(phi)*sin(psi); cos(psi)] .* vatom;

	return v_sc, Et, psi, phi
end

function trapz(x, y) # trapz performs a numerical integration after the trapezoidal rule
   if (length(y) != length(x))
       error("Vectors must be of same length")
   end
   Intresult = sum( (x[2:end] .- x[1:end-1]).*(y[2:end].+y[1:end-1]) ) / 2
   return Intresult
end
