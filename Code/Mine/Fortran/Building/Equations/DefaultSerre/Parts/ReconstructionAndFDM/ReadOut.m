% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/EnerguhtoGtou/'
Normsfile = strcat(wdir,'Norms.dat');
Energyfile = strcat(wdir,'Energy.dat');
% figure;
% Norms = importdata(Normsfile);
% 
% dxs = Norms(:,1);
% Normhs = Norms(:,2);
% NormGs = Norms(:,3);
% Normus = Norms(:,4);
% Normdus = Norms(:,5);
% 
% slope1 = 0.1*dxs.^(1);
% slope2 = 0.1*dxs.^(2);
% loglog(dxs,Normhs,'or',dxs,NormGs,'.g',dxs,Normus,'ob',dxs,Normdus,'.k' ,dxs,slope1,'-k',dxs,slope2,'-b')
% grid on
% legend('h','G', 'u','du','Slope 1', 'Slope 2','Location','northwest')

figure;
C1 = importdata(Energyfile);

dxs = C1(:,1);
C1hs = C1(:,2);
%C1Gs = C1(:,3);
%C1uhs = C1(:,4);
C1Hs = C1(:,5);


slope2 = 0.1*dxs.^(2);
slope3 = 0.1*dxs.^(3);
slope4 = 0.1*dxs.^(4);
%loglog(dxs,C1hs,'or',dxs,C1Gs,'.g',dxs,C1uhs,'ob',dxs,C1Hs, '.k',dxs,slope2,'-r',dxs,slope3,'-k',dxs,slope4,'-b')
loglog(dxs,C1hs,dxs,C1Hs, '.k',dxs,slope2,'-r',dxs,slope3,'-k',dxs,slope4,'-b')
grid on
legend('h','H','Slope 2','Slope 3', 'Slope 4','Location','northwest')



