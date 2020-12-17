% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = "./Exp/Run/NewLimrSV/10s/06/"

xhug = importdata(strcat(wdir, 'xhuGFin.dat'));
x = xhug(:,1);
h = xhug(:,2);
u = xhug(:,3);
g = xhug(:,4);

xhugA = importdata(strcat(wdir, 'CAxhuGFinA.dat'));
xA = xhugA(:,1);
hA = xhugA(:,2);
uA = xhugA(:,3);
gA = xhugA(:,4);

figure;
plot(xA,hA,'-b');
hold on;
plot(x,h,'-r');
hold off

figure;
plot(xA,gA,'-b');
hold on;
plot(x,g,'-r');
hold off

figure;
plot(xA,uA,'-b');
hold on;
plot(x,u,'-r');
hold off
