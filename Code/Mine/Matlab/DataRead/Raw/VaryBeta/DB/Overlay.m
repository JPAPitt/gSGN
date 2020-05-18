% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN/VaryBeta/ImpDispSerre/05/'
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

%figure;
%plot(x0,h0);

figure;
plot(x1,h1,'-b',x1,h1a,'--k');





