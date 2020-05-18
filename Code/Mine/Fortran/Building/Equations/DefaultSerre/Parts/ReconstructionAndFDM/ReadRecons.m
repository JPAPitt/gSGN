% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
xhGufile = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/GenInitFDMRECON/08/All.dat';

xhGu = importdata(xhGufile);

x = xhGu(:,1);
hRA = xhGu(:,2);
hR = xhGu(:,3);
GRA = xhGu(:,4);
GR = xhGu(:,5);
uRA = xhGu(:,6);
uR = xhGu(:,7);
duRA = xhGu(:,8);
duR = xhGu(:,9);

subplot(4, 3, 1);
plot(x,hRA);

subplot(4, 3, 2);
plot(x,hR);

subplot(4, 3, 3);
plot(x,hRA - hR);

subplot(4, 3, 4);
plot(x,GRA);

subplot(4, 3, 5);
plot(x,GR);

subplot(4, 3, 6);
plot(x,GRA - GR);

subplot(4, 3, 7);
plot(x,uRA);

subplot(4, 3, 8);
plot(x,uR);

subplot(4, 3, 9);
plot(x,uRA - uR);

subplot(4, 3, 10);
plot(x,duRA);

subplot(4, 3, 11);
plot(x,duR);

subplot(4, 3, 12);
plot(x,duRA - duR);
