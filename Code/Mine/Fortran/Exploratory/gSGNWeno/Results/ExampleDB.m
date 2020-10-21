% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Solver/DB21/15sCT/05/"


xhuGN = importdata(strcat(wdir, 'xhhGG.dat'));
xN = xhuGN(:,1);
hN = xhuGN(:,2);
hA = xhuGN(:,3);
GN = xhuGN(:,4);
GA = xhuGN(:,5);



figure;
plot(xN,hA);
hold on;
plot(xN,hN,'.r');
xlabel('x')
ylabel('h')


figure;
plot(xN,GA);
hold on;
plot(xN,GN,'.r');
xlabel('x')
ylabel('G')



