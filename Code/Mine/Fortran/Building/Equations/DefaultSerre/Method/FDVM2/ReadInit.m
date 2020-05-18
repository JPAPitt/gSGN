% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/';
xexfile = strcat(wdir, '08/Init.dat' );

xhGu = importdata(xexfile);

x = xhGu(:,1);
hi = xhGu(:,2);
Gi = xhGu(:,3);
ui = xhGu(:,4);


subplot(1, 3, 1);
plot(x,hi);

subplot(1, 3, 2);
plot(x,Gi);

subplot(1, 3, 3);
plot(x,ui);

