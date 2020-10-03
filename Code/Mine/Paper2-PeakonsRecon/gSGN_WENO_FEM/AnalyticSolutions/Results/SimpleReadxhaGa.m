% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Recon/Test1/00/"


xhuG = importdata(strcat(wdir, 'xhaGa.dat'));
x = xhuG(:,1);
h = xhuG(:,2);
G = xhuG(:,3);

figure;
plot(x,h);

figure;
plot(x,G);





