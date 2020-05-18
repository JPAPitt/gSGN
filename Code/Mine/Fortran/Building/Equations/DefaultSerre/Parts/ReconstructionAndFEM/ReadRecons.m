% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
whGfile = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/ReconFEM/01/ReconHG.dat';
wufile = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/ReconFEM/01/Reconu.dat';

xhG = importdata(whGfile);

xR = xhG(:,1);
hAR = xhG(:,2);
hR = xhG(:,3);
GAR = xhG(:,4);
GR = xhG(:,5);

xu = importdata(wufile);
xuR = xu(:,1);
uaR = xu(:,2);
uR = xu(:,3);

subplot(3, 3, 1);
plot(xR,hAR);

subplot(3, 3, 2);
plot(xR,hR);

subplot(3, 3, 3);
plot(xR,hAR - hR);

subplot(3, 3, 4);
plot(xR,GAR);

subplot(3, 3, 5);
plot(xR,GR);

subplot(3, 3, 6);
plot(xR,GAR - GR);

subplot(3, 3, 7);
plot(xuR,uaR);

subplot(3, 3, 8);
plot(xuR,uR);

subplot(3, 3, 9);
plot(xuR,uaR - uR);
