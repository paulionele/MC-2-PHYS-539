function Photon_Energy = energyBin(N)
%Call this function using: energyBin(N); where N is number of particles in
%simulation.

%Easy to view accuracy/distribtion with hist(Photon_Energy,24) b/c there
%are 24 energy bins in the loop. Returns an array with the exact number of
%photons at each discrete energy as specified by the energy distribution
%provided to us.

x = 0:(1E5/N):100000;
x(1) = 1; % I just didn't like the first value = 0
Photon_Energy = zeros(size(x));

for i = 1:length(x)
    
    n = x(i);
    
    if((n>0)&&(n<=2480))
        Photon_Energy(i)=0.25; % Photon Energy in MeV
    elseif((n>2480)&&(n<=15000))
        Photon_Energy(i)=0.5; % Photon Energy in MeV
    elseif((n>15000)&&(n<=27290))
        Photon_Energy(i)=0.75; % Photon Energy in MeV
    elseif((n>27290)&&(n<=37590))
        Photon_Energy(i)=1; % Photon Energy in MeV
    elseif((n>37590)&&(n<=46310))
        Photon_Energy(i)=1.25; % Photon Energy in MeV
    elseif((n>46310)&&(n<=53760))
        Photon_Energy(i)=1.5; % Photon Energy in MeV
    elseif((n>53760)&&(n<=60140))
        Photon_Energy(i)=1.75; % Photon Energy in MeV
    elseif((n>60140)&&(n<=65680))
        Photon_Energy(i)=2; % Photon Energy in MeV
    elseif((n>65680)&&(n<=70460))
        Photon_Energy(i)=2.25; % Photon Energy in MeV
    elseif((n>70460)&&(n<=74630))
        Photon_Energy(i)=2.5; % Photon Energy in MeV
    elseif((n>74630)&&(n<=78290))
        Photon_Energy(i)=2.75; % Photon Energy in MeV
    elseif((n>78290)&&(n<=81510))
        Photon_Energy(i)=3; % Photon Energy in MeV
    elseif((n>81510)&&(n<=84330))
        Photon_Energy(i)=3.25; % Photon Energy in MeV
    elseif((n>84330)&&(n<=86860))
        Photon_Energy(i)=3.5; % Photon Energy in MeV
    elseif((n>86860)&&(n<=89090))
        Photon_Energy(i)=3.75; % Photon Energy in MeV
    elseif((n>89090)&&(n<=91060))
        Photon_Energy(i)=4; % Photon Energy in MeV
    elseif((n>91060)&&(n<=92790))
        Photon_Energy(i)=4.25; % Photon Energy in MeV
    elseif((n>92790)&&(n<=94330))
        Photon_Energy(i)=4.5; % Photon Energy in MeV
    elseif((n>94330)&&(n<=95670))
        Photon_Energy(i)=4.75; % Photon Energy in MeV
    elseif((n>95670)&&(n<=96840))
        Photon_Energy(i)=5; % Photon Energy in MeV
    elseif((n>96840)&&(n<=97850))
        Photon_Energy(i)=5.25; % Photon Energy in MeV
    elseif((n>97850)&&(n<=98710))
        Photon_Energy(i)=5.5; % Photon Energy in MeV
    elseif((n>98710)&&(n<=99420))
        Photon_Energy(i)=5.75; % Photon Energy in MeV
    elseif((n>99420)&&(n<=100000))
        Photon_Energy(i)=6; % Photon Energy in MeV
    end
end


