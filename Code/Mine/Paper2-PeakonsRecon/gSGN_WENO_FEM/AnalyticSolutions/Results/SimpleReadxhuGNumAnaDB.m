% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Run/DBTest/10s/06/"


xhuGN = importdata(strcat(wdir, 'xhuGFin.dat'));
xN = xhuGN(:,1);
hN = xhuGN(:,2);
uN = xhuGN(:,3);
GN = xhuGN(:,4);

%xhuGI = importdata(strcat(wdir, 'xhuGInit.dat'));
xhuGI = importdata(strcat(wdir, 'xhuGFinA.dat'));
xI = xhuGI(:,1);
hI = xhuGI(:,2);
uI = xhuGI(:,3);
GI = xhuGI(:,4);


figure;
plot(xI,hI);
hold on;
plot(xN,hN,'.r');
xlabel('x')
ylabel('h')


figure;
plot(xI,GI);
hold on;
plot(xN,GN,'.r');
xlabel('x')
ylabel('G')

figure;
plot(xI,uI,'');
hold on;
plot(xN,uN,'.r');
xlabel('x')
ylabel('u')





