% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = './Validation/Run/SV/10s/06/';

xhug = importdata(strcat(wdir, 'xhuGInit.dat'));
x = xhug(:,1);
h = xhug(:,2);
u = xhug(:,3);
g = xhug(:,4);

xhugA = importdata(strcat(wdir, 'xhuGFinA.dat'));
xA = xhugA(:,1);
hA = xhugA(:,2);
uA = xhugA(:,3);
gA = xhugA(:,4);

figure;
plot(xA,hA,'-b');
hold on;
plot(x,h,'-r');
hold off

