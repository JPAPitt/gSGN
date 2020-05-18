% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/GenSWWE/Forced/eps1o3/';
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/GenSWWE/Forced/eps1o3/';
xexfile = strcat(wdir, '08/AllVarsEnd.dat' );

xhGu = importdata(xexfile);

x = xhGu(:,1);
hA = xhGu(:,2);
h = xhGu(:,3);
GA = xhGu(:,4);
G = xhGu(:,5);
uA = xhGu(:,6);
u = xhGu(:,7);
 

subplot(4, 3, 1);
plot(x,hA);

subplot(4, 3, 2);
plot(x,h);

subplot(4, 3, 3);
plot(x,hA - h);

subplot(4, 3, 4);
plot(x,GA);

subplot(4, 3, 5);
plot(x,G);

subplot(4, 3, 6);
plot(x,GA - G);

subplot(4, 3, 7);
plot(x,uA);

subplot(4, 3, 8);
plot(x,u);

subplot(4, 3, 9);
plot(x,uA - u);

