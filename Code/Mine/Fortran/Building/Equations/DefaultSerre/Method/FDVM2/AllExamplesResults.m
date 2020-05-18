% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir1 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta1/';
wdir0p5 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta0p5/';
wdir0 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta0/';



xhGueps1 = importdata(strcat(wdir1,'08/', 'AllVarsEnd.dat' ));
xhGueps0p5 = importdata(strcat(wdir0p5,'08/', 'AllVarsEnd.dat' ));
xhGueps0 = importdata(strcat(wdir0,'08/', 'AllVarsEnd.dat' ));

xhGueps1Init = importdata(strcat(wdir1,'08/', 'InitVars.dat' ));
xhGueps0p5Init = importdata(strcat(wdir0p5,'08/', 'InitVars.dat' ));
xhGueps0Init = importdata(strcat(wdir0,'08/', 'InitVars.dat' ));

init = xhGueps1Init;
xhGu  = xhGueps1;

x = xhGu(:,1);
hA = xhGu(:,2);
h = xhGu(:,3);
GA = xhGu(:,4);
G = xhGu(:,5);
uA = xhGu(:,6);
u = xhGu(:,7);

x0 = init(:,1);
h0 = init(:,2);
G0 = init(:,3);
u0 = init(:,4);

subplot(4, 4, 1);
plot(x0,h0);
xlabel('x (m)');
ylabel('height (m)');
title('t = 0s');


subplot(4, 4, 2);
plot(x,h);
xlabel('x (m)');
ylabel('height (m)');
title('Numerical Solution (t = 10s)');

subplot(4, 4, 3);
plot(x,hA);
xlabel('x (m)');
ylabel('height (m)');
title('Forced Solution (t = 10s)');

subplot(4, 4, 4);
plot(x,(hA - h) /hA);
xlabel('x (m)');
ylabel('h_{forced} - h_{computed}');
title('Error (t = 10s)');

