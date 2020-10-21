% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Recon/SINE/06/"


xhuGN = importdata(strcat(wdir, 'Recons.dat'));
xN = xhuGN(:,1);
hN = xhuGN(:,2);
hA = xhuGN(:,3);



figure;
plot(xN,hA);
hold on;
plot(xN,hN,'.r');
xlabel('x')
ylabel('h')






