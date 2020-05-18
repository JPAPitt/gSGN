% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/SolitonInitGCalc/';
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/MethodParts/Reconstructions/';
NormsFile = strcat(wdir, 'Norms.dat');

Out = importdata(NormsFile);

dxs = Out(:,1);
Normqs = Out(:,2);
Normdqs = Out(:,3);
Normddqs = Out(:,4);

figure;
slope1 = 0.1*dxs.^(1);
slope2 = 0.1*dxs.^(2);
loglog(dxs,Normqs,'o',dxs,Normdqs,'.',dxs,Normddqs,'o',dxs,slope1,'--k',dxs,slope2,'--b')
grid on
legend('q','dq', 'ddq','Slope 1', 'Slope 2','Location','northwest')


% Energ = importdata(EnergysFile);
% 
% dxs = Energ(:,1);
% C1h = Energ(:,2);
% C1G = Energ(:,3);
% C1uh = Energ(:,4);
% C1H = Energ(:,5);
% 
% 
% figure;
% slope1 = 0.1*dxs.^(1);
% slope2 = 0.1*dxs.^(2);
% loglog(dxs,C1h,'o',dxs,C1G,'.',dxs,C1uh,'o',dxs,C1H,'.',dxs,slope1,'--k',dxs,slope2,'--b')
% grid on
% legend('h','G', 'uh', 'H','Slope 1', 'Slope 2','Location','northwest')




