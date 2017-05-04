%MC-2 (cartersian particle position). BY ANDY OGILVY & PAUL IONELE
%Photon and basic charged particle transport simulation.

%Refer to report Problem 1 section for simulation information.

clear;clc
tic

%Called function returns arrays of energies and atten. coeffs for various.
%interactions. Note a single unique index specifies cooresponding values.
file1 = 'mass_atten.txt';
[energies,coherent,comptons,photoels,pairtrip,energytr,energyab,radifrac] = table1(file1);

mu = zeros(1,length(energies)); %total mass atten. coeff. at each energy.
for k = 1:length(energies)
    mu(k) = coherent(k) + comptons(k) + photoels(k) + pairtrip(k);
end

%Abitrary constants.
MP = 10000; %the number of primary photons to simulate.
NP = 7;  %this defines the number of properties for each particle.
ec = 0.01; %photon and electron cut-off energy in MeV.
ctheta = 0:0.01:pi; %range of angles between [0,180] for Compton.

%Pre-allocating a large zero-matrix that will fill with the values of dose
%and position corresponding to either photons or electrons/positrons.
blockSizeP = MP; %theoretical max number of photons going below cutoff.
blockSizeE = 150*MP; %roughly arbirarily set size.
blockSizeET = 100*MP; %enough to contain all possible photon energy transfer events

%'positionDosePhoton' and 'positionDoseElectron' are the energy desposition
%arrays. They contain a record of the location of each energy deposition
%event within the tank volume. The deposition events for photons come from
%energies below cut-off, for electrons/positrons they come from their CSDA
%nature. [x,y,z,energy]
positionDosePhoton = zeros(blockSizeP,4); %photons
positionDoseElectron = zeros(blockSizeE,4); %electrons and positrons
positionEnergyTransfer = zeros(blockSizeET,4); %energy transfer


%Counters to ensure the above dose-dep arrays fill sequentially.
count1 = 1;
count2 = 1;
count3 = 1;

%Setting electron step-size and energy lost per step.
elecStepSize = 0.01; %elecStepSize = 0.01cm, or 0.1mm
energyPerStep = elecStepSize*2; %energy deposited per step

%%%PARTICLE STACK. MP dim is undetermined but is at least MP (initial #).
%Ideally we would like to preallocate as much space as possible, but this
%stack will grow as new particles are created. New particle arrays are
%concatenated to this array. This is likely the performance limiter.
%%%%%Current stack configuration: particles #s in x, properties in y
%[type,energy,x,y,z,theta,phiV]
% Where x,y,z are cartesian coordinates for particle position.
% Where theta, phiV are spherical unit vectors for particle direction.
stack = zeros(length(MP),length(NP));

%Energy spectrum of incoming photons (1 x length(MP)) array.
spectrum = energyBin(MP);

%Determining random position on water surface, working in cylindrical coor.
r1 = (-2*log(rand(MP,1))).^.5; %r coordinate (non uniform; Gaussian)
r2 = rand(MP,1)*2*pi; %phi coordinate (uniform random on [0,2*pi])
for i = 1:length(r1)
    %Loop is to reject any photons that are greater than 5cm from axis.
    if r1(i) > 5
        while r1(i) > 5
            %Keep trying for random number until satisfied. Replace
            %rejected value so we still start with MP photons.
            r1(i) = (-2*log(rand(1))).^.5;
        end
    end
end

%Construct initial stack given array of MP photons at defined energies.
%Photon energies set by distribution, positions are set random at surface.
for i = 1:MP
    %Set particle type.
    stack(i,1) = 0; %type 0 = photon, 1 = electron, 2 = positron
    
    %Get and set photon energy.
    stack(i,2) = spectrum(i);
    
    %Determine position on water surface randomly.
    stack(i,3) = r1(i); %rho
    stack(i,4) = r2(i); %phiC, random on (0,360] (angle in x-y plane)
    stack(i,5) = 0; %z starts at zero (surface)
    
    %Initial photon directions are all down!!! No divergence from source.
    %This is one of the approximations we made in this simulation.
    stack(i,6) = 0; %theta, random 0 to 180 (angle measured from (-)z-axis)
    stack(i,7) = 0; %phiV, random 0 to 360 (angle measured in xy-plane)
    
end

%Particle positions from now on will be in cartesian coordinates.
for i = 1:MP
    %Re-express particle positions on the water surface in cartesian.
    stack(i,3) = r1(i)*cos(r2(i)); %coordinate x
    stack(i,4) = r1(i)*sin(r2(i)); %coordinate y
end
dispCount = 1;
while size(stack,1) ~= 0
    %While the stack is not empty...
    for i = 1:1 %vestigial loop, literally the tailbone or appendix of loops
        %For every particle in the stack...
        if stack(i,1) == 0
            %%%START PHOTONS --- Bit value 0 for photon.
            %Sample interaction type.
            if stack(i,2) >= ec %while photon energy is greater than cutoff energy
                
                %Linear interpolation between linear attenuation coefficients.
                %DO NOT REUSE r3!!!
                mu_interp = interp1(energies, mu, stack(i,2));
                r3 = (-1/mu_interp)*log(rand(1)); %randomly sampled distance.
                
                %Location of interaction after traveling mean free path with
                %angles phiV and theta. The new position of the particle is
                %determined based on the the trajectory of the previous and the
                %random distance travelled computed above. The sin and cos are
                %due to conversion from spherical to cartesian coordinates.
                %[type,energy,x,y,z,theta,phiV] <--- reference
                stack(i,3) = stack(i,3) + r3*sin(stack(i,6))*cos(stack(i,7));
                stack(i,4) = stack(i,4) + r3*sin(stack(i,6))*sin(stack(i,7));
                stack(i,5) = stack(i,5) + r3*cos(stack(i,6));
                
                %Check if current trajectory and new calculated distance take
                %photon outside of volume of interest. If yes, remove photon
                %from stack. If no, continue to determine if energy above
                %threshold cutoff.
                if ((abs(stack(i,3)) > 15) || (abs(stack(i,4)) > 15))...
                        || ((stack(i,5) > 50) || (stack(i,5) <= 0))
                    stack(i,:) = []; %Clear stack of current photon.
                    %disp('photon removed, outside tank')
                    continue %continue on to the next particle in the stack.
                end
                
                %Computing interaction probabilities; use mu_interp from
                %above. Normalize individual interaction probabilites by
                %mu_interp to obtain relative probabilities.
                %Linear interpolation between mass attenuation coefficients.
                coherent_interp = interp1(energies, coherent, stack(i,2));
                comptons_interp = interp1(energies, comptons, stack(i,2));
                photoels_interp = interp1(energies, photoels, stack(i,2));
                pairtrip_interp = interp1(energies, pairtrip, stack(i,2));
                
                %Probabilities for interactions, used to select interaction
                p1 =      coherent_interp/mu_interp;
                p2 = p1 + comptons_interp/mu_interp;
                p3 = p2 + photoels_interp/mu_interp;
                p4 = p3 + pairtrip_interp/mu_interp;
                
                r1 = rand(1); %generate uniform random deviate on [0,1].
                
                if (r1 < p1)
                    %COHERENT SCATTER.
                    
                    %Assuming uniform randomly distributed new photon
                    %direction and update information in the particle's
                    %information.
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack(i,6) = rand(1)*pi; %theta, random 0 to 180 (angle measured from (-)z-axis)
                    stack(i,7) = rand(1)*2*pi; %phiV, random 0 to 360 (angle measured on x-y plane)
                    
                    positionEnergyTransfer(count3,:) = [stack(i,3), stack(i,4), stack(i,5), 0];
                    
                    %Count3 just makes sure positionDosePhoton fills
                    %sequentially without skipping rows
                    count3 = count3 + 1;
                    
                elseif (r1 >= p1) && (r1 < p2)
                    %COMPTON SCATTER.
                    
                    %The following is to determine new photon energy.
                    %First, an array of photon energies ep, based on photon
                    %scatter angles [0,180] with fine sampling; theta
                    %defined in header. stack(i,2) is photon current energy
                    ep = stack(i,2)./(1+(stack(i,2)/0.51099)*(1-cos(ctheta)));
                    
                    %Here, we build the differential compton cross section
                    %for the e MeV photon energy. This is again an array. We
                    %will sample from this array.
                    d_cr = (2*pi*sin(ctheta)).*((2.8179*10^-15)^2)/2.*...
                        ((ep/stack(i,2)).^2).*(stack(i,2)./ep + ep/stack(i,2)...
                        - (sin(ctheta).^2));
                    d_cr_max = max(d_cr); %returns maximum value of d_cr; no index
                    
                    while 1
                        %Using the 'rejection method' to determine energy.
                        %This loop will continue as long as r4 is rejected.
                        %Loop terminates by 'break'.
                        r4 = rand(1)*pi; %uniform random deviate on [0,pi]
                        r5 = rand(1)*d_cr_max; %[0,max_cross_section]
                        
                        %Get scattered photon energy ep2 at the random
                        %angle r2. Evaluate d_cr again but for the random
                        %angle and cooresponding photon energy ep2 ->
                        %d_cr2.
                        ep2 = stack(i,2)/(1+(stack(i,2)/0.51099)*(1-cos(r4)));
                        d_cr2 = (2*pi*sin(r4))*((2.8179*10^-15)^2)/2*...
                            ((ep2/stack(i,2))^2)*(stack(i,2)/ep2 + ep2/stack(i,2) - (sin(r4)^2));
                        
                        if r5 <= d_cr2
                            %Accept r5 random angle, continue on to
                            %calculate the kinematic parameters of the
                            %scattered photon and Compton electron.
                            %Determine angles.
                            
                            %Electron kinetic energy = initial photon e - scattered e
                            TE = stack(i,2) - ep2;
                            
                            positionEnergyTransfer(count3,:) = [stack(i,3), stack(i,4), stack(i,5), TE];
                            
                            %Count3 just makes sure positionDosePhoton fills
                            %sequentially without skipping rows
                            count3 = count3 + 1;
                            
                            %Photon changes.
                            original_photon_energy = stack(i,2); %still need to old energy for electron energry calculation
                            stack(i,2) = ep2; %update photon energy to new scattered photon energy
                            rotAng = rand(1)*2*pi; %eta, rotational angle around initial photon trajectory, assuming uniform distn
                            photAng = r4;
                            photDir = directionCosine(stack(i,6),stack(i,7),rotAng,photAng); %returns the correct photon direction
                            stack(i,6) = photDir(1); %update photon direction
                            stack(i,7) = photDir(2);
                            
                            %Electron changes.
                            elecAng = atan(cot(photAng/2)/(1+original_photon_energy/0.511));
                            elecDir = directionCosine(stack(i,6),stack(i,7),rotAng,elecAng); %correct electron direction
                            
                            %Add new electron to stack.
                            stack = [stack; [1, TE, stack(i,3),stack(i,4),stack(i,5),elecDir(1),elecDir(2)]];
                            
                            break %breaks out of while loop
                        else
                            %Reject r2 random angle, run while loop again.
                        end
                    end
                elseif (r1 >= p2) && (r1 < p3)
                    %PHOTOELECTRIC.
                    
                    %Energy transferred to electron.
                    TE = stack(i,2) - 543.1E-6; %binding energy = 543eV
                    
                    %Create an electron with same direction as photon.
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack = [stack; [1, TE, stack(i,3),stack(i,4),stack(i,5),stack(i,6),stack(i,7)]];
                    
                    positionEnergyTransfer(count3,:) = [stack(i,3), stack(i,4), stack(i,5), TE];
                    
                    %Count3 just makes sure positionDosePhoton fills
                    %sequentially without skipping rows
                    count3 = count3 + 1;
                    
                    stack(i,:) = []; %delete photon from stack
                    %disp('photon removed by PE')
                    
                else
                    %PAIR + TRIPLET PRODUCTION.
                    
                    %Sample angles from gaussian distribution with mean
                    %value at original photon angle/trajectory. Energy
                    %transferred to photons is:
                    TE = stack(i,2) - 1.022;
                    
                    meanAng = 0.511/(TE/2);
                    elecAng = acos(1-2*rand(1))- pi + meanAng; %Chi, guassian distn mean = pi, subtract pi and add the mean we want?
                    posAng = -acos(1-2*rand(1))- pi + meanAng; %Chi, guassian distn mean = pi, subtract pi and add the mean we want?
                    rotAng = rand(1)*2*pi; %Eta, rotational angle, assuming uniform distribution
                    
                    %Returning the correct directions for photons & elecs.
                    posDir = directionCosine(stack(i,6),stack(i,7),rotAng,posAng);
                    elecDir = directionCosine(stack(i,6),stack(i,7),rotAng,elecAng);
                    
                    %Popping a positron and electron onto the stack.
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack = [stack; [2, TE/2, stack(i,3),stack(i,4),stack(i,5),posDir(1),posDir(2)]];
                    stack = [stack; [1, TE/2, stack(i,3),stack(i,4),stack(i,5),elecDir(1),elecDir(2)]];
                    
                    positionEnergyTransfer(count3,:) = [stack(i,3), stack(i,4), stack(i,5), TE];
                    
                    %Count3 just makes sure positionDosePhoton fills
                    %sequentially without skipping rows
                    count3 = count3 + 1;
                    
                    stack(i,:) = []; %delete photon from stack
                    %disp('photon removed by PP')
                end
                %At this point, photon may have new energy dependent on the
                %interaction that it underwent, it may have a new
                %direction, it may have been totally absorbed and or new
                %particles may have been produced.
            else
                %Photon energy is below cutoff energy. Deposit energy
                %remaining at current location. Delete photon from stack.
                
                %[type,energy,x,y,z,theta,phiV] <--- reference
                positionDosePhoton(count1,:) = [stack(i,3), stack(i,4), stack(i,5), stack(i,2)];
                
                %Count1 just makes sure positionDosePhoton fills
                %sequentially without skipping rows
                count1 = count1 + 1;
                
                stack(i,:) = []; %delete photon from stack
                %disp('photon removed, low energy')
            end
            %%%END PHOTONS
        else
            %%%START PARTICLES
            %Non-zero bit for charged particle. Alg. for charged particles.
            if stack(i,1) == 1
                %%%START ELECTRONS --- Bit value 1 for electron.
                
                %Check if current trajectory and new calculated distance
                %take electron outside of tank volume. If yes, remove
                %electron from stack. If no, continue to determine if
                %energy above threshold cutoff.
                if ((abs(stack(i,3)) > 15) || (abs(stack(i,4)) > 15))...
                        || ((stack(i,5) > 50) || (stack(i,5) <= 0))
                    stack(i,:) = []; %delete electron from stack
                    %disp('electron removed, outside tank')
                    continue %Continue on to the next particle in the stack.
                end
                
                if stack(i,2) > energyPerStep
                    %If the energy of the electron is greater than the
                    %energy lost per step, update energy and position of
                    %electron. %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack(i,3) = stack(i,3) + energyPerStep*sin(stack(i,6))*cos(stack(i,7));
                    stack(i,4) = stack(i,4) + energyPerStep*sin(stack(i,6))*sin(stack(i,7));
                    stack(i,5) = stack(i,5) + energyPerStep*cos(stack(i,6));
                    
                    %[x,y,z,energy] <--reference
                    positionDoseElectron(count2,:) = [stack(i,3), stack(i,4), stack(i,5), energyPerStep];
                    
                    %Count2 just makes sure positionDoseElectron fills
                    %sequentially without skipping rows.
                    count2 = count2 + 1;
                    
                    %Update electron energy.
                    stack(i,2) = stack(i,2) - energyPerStep;
                    
                else
                    %Electron energy remaining is less than energy lost in
                    %a step. Deposit remaining at the final location and
                    %calculate total distance electron travels along
                    %trajectory. The new distance will be less than the
                    %previous travelled distances.
                    
                    %Distance in cm, determined from given rate of energy
                    %loss: 2MeV/cm.
                    elecDistance = stack(i,2)/2;
                    
                    %New electron position based on elecDistance step size.
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack(i,3) = stack(i,3) + elecDistance*sin(stack(i,6))*cos(stack(i,7));
                    stack(i,4) = stack(i,4) + elecDistance*sin(stack(i,6))*sin(stack(i,7));
                    stack(i,5) = stack(i,5) + elecDistance*cos(stack(i,6));
                    
                    %[x,y,z,energy] <--reference
                    positionDoseElectron(count2,:) = [stack(i,3), stack(i,4), stack(i,5), stack(i,2)];
                    
                    %Count2 just makes sure positionDoseElectron fills
                    %sequentially without skipping rows.
                    count2 = count2 + 1;
                    
                    stack(i,:) = []; %delete electron from stack
                    %disp('electron removed, energy depleted')
                    
                end
                %%% END ELECTRONS.
            else
                %%% START POSITRONS --- Bit value 2 for positron.
                
                %Check if current trajectory and new calculated distance
                %take electron outside of tank volume.
                if ((abs(stack(i,3)) > 15) || (abs(stack(i,4)) > 15))...
                        || ((stack(i,5) > 50) || (stack(i,5) <= 0))
                    stack(i,:) = []; %Clear stack of current positron.
                    %disp('positron removed, outside tank')
                    continue %continue on to the next particle in the stack
                end
                
                if stack(i,2) > energyPerStep
                    %If the energy of the positron is greater than the
                    %energy lost per step, update energy and position of
                    %positron. Same step energy lost for electrons and
                    %positrons!
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack(i,3) = stack(i,3) + energyPerStep*sin(stack(i,6))*cos(stack(i,7));
                    stack(i,4) = stack(i,4) + energyPerStep*sin(stack(i,6))*sin(stack(i,7));
                    stack(i,5) = stack(i,5) + energyPerStep*cos(stack(i,6));
                    
                    %[x,y,z,energy] <--reference
                    positionDoseElectron(count2,:) = [stack(i,3), stack(i,4), stack(i,5), energyPerStep];
                    
                    %Count2 just makes sure positionDoseElectron fills
                    %sequentially without skipping rows.
                    count2 = count2 + 1;
                    
                    %Update positron energy.
                    stack(i,2) = stack(i,2) - energyPerStep;
                    
                else
                    %Positron energy remaining is less than energy lost in
                    %a step. Deposit remaining at the final location and
                    %calculate total distance positron travels along
                    %trajectory. The new distance will be less than the
                    %previous travelled distances.
                    
                    %Distance in cm, determined from given rate of energy
                    %loss: 2MeV/cm.
                    posDistance = stack(i,2)/2;
                    
                    %New electron position based on posDistance step size.
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack(i,3) = stack(i,3) + posDistance*sin(stack(i,6))*cos(stack(i,7));
                    stack(i,4) = stack(i,4) + posDistance*sin(stack(i,6))*sin(stack(i,7));
                    stack(i,5) = stack(i,5) + posDistance*cos(stack(i,6));
                    
                    %[x,y,z,energy] <--reference
                    positionDoseElectron(count2,:) = [stack(i,3), stack(i,4), stack(i,5), stack(i,2)];
                    
                    %Count2 just makes sure positionDoseElectron fills
                    %sequentially without skipping rows.
                    count2 = count2 + 1;
                    
                    %Create two photons with opposite direction, with
                    %'common' random rotational angle.
                    rotAng = rand(1)*2*pi;
                    %One photon is created at 90 degrees to the positron
                    %direction. The other is opposite to the first at 270
                    %degrees.
                    photAng1 = pi/2;
                    photAng2 = 3*pi/2;
                    
                    %Getting the correct photon directions.
                    photDir1 = directionCosine(stack(i,6),stack(i,7),rotAng,photAng1);
                    photDir2 = directionCosine(stack(i,6),stack(i,7),rotAng,photAng2);
                    
                    %Add photon 1 to stack.
                    %Photon energies are well-defined: 0.511 MeV
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack(i,6) = photDir1(1);
                    stack(i,7) = photDir1(2);
                    stack = [stack; [0, 0.511, stack(i,3),stack(i,4),stack(i,5),stack(i,6),stack(i,7)]];
                    
                    %Add photon 2 to stack.
                    %[type,energy,x,y,z,theta,phiV] <--- reference
                    stack(i,6) = photDir2(1);
                    stack(i,7) = photDir2(2);
                    stack = [stack; [0, 0.511, stack(i,3),stack(i,4),stack(i,5),stack(i,6),stack(i,7)]];
                    
                    stack(i,:) = []; %delete positron from stack
                    %disp('positron removed, annihilated')
                end
                %%% END POSITRONS
            end
            %%% END PARTICLES
        end
        
        dispCount = dispCount + 1;
        if (rem(dispCount,10000)==0)
            disp(size(stack,1));
        end
        
    end
end
disp('stack finished');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Binning dose deposition events into voxels...

%Remove excess preallocated zeros in positionDoseElectron array.
positionDoseElectron(all(positionDoseElectron==0,2),:) = [];
positionEnergyTransfer(all(positionEnergyTransfer==0,2),:) = [];

%Change deposition locations to 0to30cm x 0to30cm, from -15to15cm x -15to15cm
positionDoseElectron(:,1) = positionDoseElectron(:,1)+15;
positionDoseElectron(:,2) = positionDoseElectron(:,2)+15;
positionEnergyTransfer(:,1) = positionEnergyTransfer(:,1)+15;
positionEnergyTransfer(:,2) = positionEnergyTransfer(:,2)+15;

%Some electrons may have one track point outside the tank (reason unknown),
%remove those points from the dose deposition array.
positionDoseElectron(any(positionDoseElectron < 0,2),:) = [];
positionEnergyTransfer(any(positionEnergyTransfer < 0,2),:) =[];

% figure(1)
% scatter(positionDoseElectron(:,1),positionDoseElectron(:,2))

%Voxel-based grids to store dose accumumated in each voxel and total number
%of interactions within that voxel. May add grids to this as required.
grid_dose = zeros(61,61,101);
grid_inte = zeros(61,61,101);

%Multiply x,y,z position from 1 cm^3 to 0.5 cm^3.
positionDoseElectron2 = [positionDoseElectron(:,1:3)*2,positionDoseElectron(:,4)];
positionEnergyTransfer2 = [positionEnergyTransfer(:,1:3)*2,positionEnergyTransfer(:,4)];

%Round down to nearest integer to get voxel location for deposition event.
voxels_dose = [floor(positionDoseElectron2(:,1:3))+1,positionDoseElectron(:,4)];


for i = 1:length(voxels_dose)
    x = voxels_dose(i,1);
    y = voxels_dose(i,2);
    z = voxels_dose(i,3);
    
    grid_dose(y,x,z) = grid_dose(y,x,z) + voxels_dose(i,4);
end

voxels_ET = [floor(positionEnergyTransfer2(:,1:3))+1,positionEnergyTransfer(:,4)];
for i = 1:length(voxels_ET)
    x = voxels_ET(i,1);
    y = voxels_ET(i,2);
    z = voxels_ET(i,3);
    
    grid_inte(y,x,z) = grid_inte(y,x,z) + voxels_ET(i,4);
end

grid_dose_prime = grid_dose;
grid_inte_prime = grid_inte;
    
depth = 1:0.5:51;
x = 1:0.5:31;
caxDose = grid_dose_prime(30,30,:) + grid_dose_prime(31,30,:) + grid_dose_prime(30,31,:) + grid_dose_prime(31,31,:);
caxKerma = grid_inte_prime(30,30,:) + grid_inte_prime(31,30,:) + grid_inte_prime(30,31,:) + grid_inte_prime(31,31,:);

caxDose = squeeze(caxDose);
caxKerma = squeeze(caxKerma);

figure(1)
plot(depth,caxDose)
xlabel('Depth(cm)')
ylabel('Dose(Arbitrary)')
title('Depth Dose Along the Z-axis')
grid on

xdepth1dot5cm = grid_dose_prime(:,30,3) + grid_dose_prime(:,31,3) + grid_dose_prime(:,30,4) + grid_dose_prime(:,31,4);
ydepth1dot5cm = grid_dose_prime(30,:,3) + grid_dose_prime(31,:,3) + grid_dose_prime(30,:,4) + grid_dose_prime(31,:,4);

xdepth5cm = grid_dose_prime(:,30,10) + grid_dose_prime(:,31,10) + grid_dose_prime(:,30,11) + grid_dose_prime(:,31,11);
ydepth5cm = grid_dose_prime(30,:,10) + grid_dose_prime(31,:,10) + grid_dose_prime(30,:,11) + grid_dose_prime(31,:,11);

xdepth10cm = grid_dose_prime(:,30,20) + grid_dose_prime(:,31,20) + grid_dose_prime(:,30,21) + grid_dose_prime(:,31,21);
ydepth10cm = grid_dose_prime(30,:,20) + grid_dose_prime(31,:,20) + grid_dose_prime(30,:,21) + grid_dose_prime(31,:,21);

xdepth20cm = grid_dose_prime(:,30,40) + grid_dose_prime(:,31,40) + grid_dose_prime(:,30,41) + grid_dose_prime(:,31,41);
ydepth20cm = grid_dose_prime(30,:,40) + grid_dose_prime(31,:,40) + grid_dose_prime(30,:,41) + grid_dose_prime(31,:,41);

figure(2)
plot(x,xdepth1dot5cm,'DisplayName','z = 1.5cm')
hold on
plot(x,xdepth5cm, 'DisplayName','z = 5cm')
hold on
plot(x,xdepth10cm,'DisplayName','z = 10cm')
hold on 
plot(x,xdepth20cm,'DisplayName','z = 20cm')
hold off
xlabel('x-Distance(cm)')
ylabel('Dose(Arbitrary)')
title('Depth Dose Along the X-axis at Various Depths')
legend('show')
grid on

figure(3)
plot(x,ydepth1dot5cm,'DisplayName','z = 1.5cm')
hold on
plot(x,ydepth5cm,'DisplayName','z = 5cm')
hold on
plot(x,ydepth10cm,'DisplayName','z = 10cm')
hold on 
plot(x,ydepth20cm,'DisplayName','z = 20cm')
hold off
xlabel('y-Distance(cm)')
ylabel('Dose(Arbitrary)')
title('Depth Dose Along the Y-axis at Various Depths')
legend('show')
grid on

figure(5)
plot(depth,mat2gray(caxDose),'DisplayName','Depth Dose')
hold on
plot(depth,mat2gray(caxKerma),'DisplayName','Depth Kerma')
hold off
legend('show')
xlabel('Depth(cm)')
ylabel('Dose(Arbitrary)')
title('Normalized Depth Dose , and Normalized Depth Kerma Along the Z-axis')
grid on

toc