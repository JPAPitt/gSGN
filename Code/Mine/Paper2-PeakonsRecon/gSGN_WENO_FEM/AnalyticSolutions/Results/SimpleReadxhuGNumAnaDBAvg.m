% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Run/DBTests/5s/06/"


xhuGN = importdata(strcat(wdir, 'xhGAvg.dat'));
xA = xhuGN(:,1);
hA = xhuGN(:,2);
GA = xhuGN(:,3);

%xhuGI = importdata(strcat(wdir, 'xhuGInit.dat'));
xhuGI = importdata(strcat(wdir, 'xhuGFinA.dat'));
xI = xhuGI(:,1);
hI = xhuGI(:,2);
uI = xhuGI(:,3);
GI = xhuGI(:,4);


figure;
plot(xI,hI);
hold on;
plot(xA,hA,'.r');
xlabel('x')
ylabel('h')


figure;
plot(xI,GI);
hold on;
plot(xA,GA,'.r');
xlabel('x')
ylabel('G')






