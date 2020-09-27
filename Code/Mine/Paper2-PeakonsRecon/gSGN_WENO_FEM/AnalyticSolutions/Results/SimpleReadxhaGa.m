% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Recon/Test1/05/"


xhuG = importdata(strcat(wdir, 'xhuGInit.dat'));
x = xhuG(:,1);
h = xhuG(:,2);
u = xhuG(:,3);
G = xhuG(:,4);

figure;
plot(x,h);

figure;
plot(x,u);

figure;
plot(x,G);




