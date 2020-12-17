% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/AnalyticSolutions/SerreSoliton/07"


xhuGN = importdata(strcat(wdir, '1.dat'));
xN = xhuGN(:,1);
hN = xhuGN(:,2);
uN = xhuGN(:,3);



figure;
plot(xN,hN);
xlabel('x')
ylabel('h')




