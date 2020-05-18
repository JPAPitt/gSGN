% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/SolitonInitGCalc/';
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/SolitonInitGAna/';
NormsFile = strcat(wdir, 'Norms.dat');
EnergysFile = strcat(wdir, 'Energy.dat');

Out = importdata(NormsFile);

dxs = Out(:,1);
Normhs = Out(:,2);
NormGs = Out(:,3);
Normus = Out(:,4);

figure;
slope1 = 0.1*dxs.^(1);
slope2 = 0.1*dxs.^(2);
loglog(dxs,Normhs,'or',dxs,NormGs,'.g',dxs,Normus,'ob',dxs,slope1,'-k',dxs,slope2,'-b')
grid on
legend('h','G', 'u','Slope 1', 'Slope 2','Location','northwest')


Energ = importdata(EnergysFile);

dxs = Energ(:,1);
C1h = Energ(:,2);
C1G = Energ(:,3);
C1uh = Energ(:,4);
C1H = Energ(:,5);


figure;
slope1 = 0.1*dxs.^(1);
slope2 = 0.1*dxs.^(2);
loglog(dxs,C1h,'or',dxs,C1G,'.g',dxs,C1uh,'ob',dxs,C1H,'.k',dxs,slope1,'-k',dxs,slope2,'-b')
grid on
legend('h','G', 'uh', 'H','Slope 1', 'Slope 2','Location','northwest')




