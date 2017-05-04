function [energies,coherent,comptons,photoels,pairtrip,energytr,energyab,radifrac] = table1(file1)
%Function to call when discrete energies are desired.
%Note: remove semi-colon at the end of each line.

%clear all;clc
fid = fopen(file1,'r');

%contents = ['Energy (MeV)' 'Coherent (cm^2g^-1)' 'Compton (cm^2g^-1)' 
%'Photoelectric (cm^2g^-1)' 'Pair+Triplet (cm^2g^-1)' 
%'Energy-Transfer (cm^2g^-1)' 'Energy-Absorption(cm^2g^-1)' '1-g'];

energies = [];
coherent = [];
comptons = [];
photoels = [];
pairtrip = [];
energytr = [];
energyab = [];
radifrac = [];

tchar = fgetl(fid);
tsplt = strsplit(tchar,' '); %strsplit returns cell type, with sub char type

st1 = str2double(tsplt(1));
st2 = str2double(tsplt(2));
st3 = str2double(tsplt(3));
st4 = str2double(tsplt(4));
st5 = str2double(tsplt(5));
st6 = str2double(tsplt(6));
st7 = str2double(tsplt(7));
st8 = str2double(tsplt(8));

energies = [energies,st1];
coherent = [coherent,st2];
comptons = [comptons,st3];
photoels = [photoels,st4];
pairtrip = [pairtrip,st5];
energytr = [energytr,st6];
energyab = [energyab,st7];
radifrac = [radifrac,st8];

while ischar(tchar)
    tchar = fgetl(fid);
    if ischar(tchar) == 1
        tsplt = strsplit(tchar,' '); %strsplit returns cell type, with sub char type
        
        st1 = str2double(tsplt(1));
        st2 = str2double(tsplt(2));
        st3 = str2double(tsplt(3));
        st4 = str2double(tsplt(4));
        st5 = str2double(tsplt(5));
        st6 = str2double(tsplt(6));
        st7 = str2double(tsplt(7));
        st8 = str2double(tsplt(8));

        energies = [energies,st1];
        coherent = [coherent,st2];
        comptons = [comptons,st3];
        photoels = [photoels,st4];
        pairtrip = [pairtrip,st5];
        energytr = [energytr,st6];
        energyab = [energyab,st7];
        radifrac = [radifrac,st8];
    end
end

fclose(fid);


