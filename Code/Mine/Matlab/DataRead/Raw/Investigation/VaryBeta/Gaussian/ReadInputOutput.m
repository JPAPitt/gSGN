% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

%gSGNForcedLimhG
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/RegSWWE/Gaussian/timeseries/06/';

initfile = strcat(wdir, 'InitVal.dat' );
endfile = strcat(wdir, 'EndVals.dat' );

xhGuinit = importdata(initfile);
xhGuend = importdata(endfile);

x0 = xhGuinit(:,1);
h0 = xhGuinit(:,2);
G0 = xhGuinit(:,3);
u0 = xhGuinit(:,4);

x1 = xhGuend(:,1);
h1 = xhGuend(:,2);
G1 = xhGuend(:,3);
u1 = xhGuend(:,4);



subplot(3, 4, 1);
plot(x0,h0);

subplot(3, 4, 2);
plot(x1,h1);

subplot(3, 4, 5);
plot(x0,G0);

subplot(3, 4, 6);
plot(x1,G1);

subplot(3, 4, 9);
plot(x0,u0);

subplot(3, 4, 10);
plot(x1,u1);