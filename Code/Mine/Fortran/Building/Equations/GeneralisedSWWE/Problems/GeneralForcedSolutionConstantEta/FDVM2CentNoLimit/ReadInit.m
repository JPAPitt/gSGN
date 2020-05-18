% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/GenSWWE/Forced/eps1o3/';
xexfile = strcat(wdir, '08/InitVars.dat' );

xhGu = importdata(xexfile);

x = xhGu(:,1);
h = xhGu(:,2);
G = xhGu(:,3);
u = xhGu(:,4);
 

subplot(1, 3, 1);
plot(x,h);

subplot(1, 3, 2);
plot(x,G);

subplot(1, 3, 3);
plot(x,u);


