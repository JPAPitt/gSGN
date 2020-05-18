% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/MethodParts/Reconstructions/';
xexfile = strcat(wdir, '10/All.dat' );

xhGu = importdata(xexfile);

x = xhGu(:,1);
qA = xhGu(:,2);
q = xhGu(:,3);
dqA = xhGu(:,4);
dq = xhGu(:,5);
ddqA = xhGu(:,6);
ddq = xhGu(:,7);
 

subplot(4, 3, 1);
plot(x,qA);

subplot(4, 3, 2);
plot(x,q);

subplot(4, 3, 3);
plot(x,qA - q);

subplot(4, 3, 4);
plot(x,dqA);

subplot(4, 3, 5);
plot(x,dq);

subplot(4, 3, 6);
plot(x,dqA - dq);

subplot(4, 3, 7);
plot(x,ddqA);

subplot(4, 3, 8);
plot(x,ddq);

subplot(4, 3, 9);
plot(x,ddqA - ddq);

dderr = ddqA - ddq;

