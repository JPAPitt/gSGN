% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/SWWE/Dambreak/'
initfile = strcat(wdir, 'InitVal.dat' );
endfile = strcat(wdir, 'EndVals.dat' );
endanafile = strcat(wdir, 'EndAnaVals.dat' );

xhGuinit = importdata(initfile);
xhGuend = importdata(endfile);
xhGuendana = importdata(endanafile);

x0 = xhGuinit(:,1);
h0 = xhGuinit(:,2);
G0 = xhGuinit(:,3);
u0 = xhGuinit(:,4);

x1 = xhGuend(:,1);
h1 = xhGuend(:,2);
G1 = xhGuend(:,3);
u1 = xhGuend(:,4);

x1a = xhGuendana(:,1);
h1a = xhGuendana(:,2);
G1a = xhGuendana(:,3);
u1a = xhGuendana(:,4);


subplot(3, 4, 1);
plot(x0,h0);

subplot(3, 4, 2);
plot(x1,h1);

subplot(3, 4, 3);
plot(x1,h1a);

subplot(3, 4, 4);
plot(x1,(h1 - h1a));

subplot(3, 4, 5);
plot(x0,G0);

subplot(3, 4, 6);
plot(x1,G1);

subplot(3, 4, 7);
plot(x1,G1a);

subplot(3, 4, 8);
plot(x1,(G1 - G1a));

subplot(3, 4, 9);
plot(x0,u0);

subplot(3, 4, 10);
plot(x1,u1);

subplot(3, 4, 11);
plot(x1,u1a);

subplot(3, 4, 12);
plot(x1,(u1 - u1a));




