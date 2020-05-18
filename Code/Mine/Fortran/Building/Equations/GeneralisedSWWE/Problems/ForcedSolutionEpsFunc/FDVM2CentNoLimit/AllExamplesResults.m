% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/GenSWWE/ForcedEpsFunc/expTimeDeriv1/';
exnum = '08/';

xhGu = importdata(strcat(wdir,exnum , 'AllVarsEnd.dat' ));
init = importdata(strcat(wdir,exnum , 'InitVars.dat' ));

x = xhGu(:,1);
hA = xhGu(:,2);
h = xhGu(:,3);
GA = xhGu(:,4);
G = xhGu(:,5);
uA = xhGu(:,6);
u = xhGu(:,7);
e = xhGu(:,8);

x0 = init(:,1);
h0 = init(:,2);
G0 = init(:,3);
u0 = init(:,4);
e0 = init(:,5);

figure;
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
plot(x,((hA - h) ./ hA));
xlabel('x (m)');
ylabel(' ^{(h_{forced} - h_{computed})}/_{h_{forced}}');
title('Error (t = 10s)');

subplot(4, 4, 5);
plot(x0,u0);
xlabel('x (m)');
ylabel('u (m/s)');
title('t = 0s');

subplot(4, 4, 6);
plot(x,u);
xlabel('x (m)');
ylabel('u (m/s)');
title('Numerical Solution (t = 10s)');

subplot(4, 4,7);
plot(x,uA);
xlabel('x (m)');
ylabel('u (m/s)');
title('Forced Solution (t = 10s)');

subplot(4, 4, 8);
plot(x,(uA - u));
xlabel('x (m)');
ylabel(' ^{(u_{forced} - u_{computed})}');
title('Error (t = 10s)');

subplot(4, 4, 9);
plot(x0,G0);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('t = 0s');

subplot(4, 4, 10);
plot(x,G);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Numerical Solution (t = 10s)');

subplot(4, 4,11);
plot(x,GA);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Forced Solution (t = 10s)');

subplot(4, 4, 12);
plot(x,(GA - G));
xlabel('x (m)');
ylabel(' ^{(G_{forced} - G_{computed})}');
title('Error (t = 10s)');

subplot(4, 4, 13);
plot(x0,e0);
xlabel('x (m)');
ylabel('\epsilon');
title('Initial Profile (t = 0s)');

subplot(4, 4, 14);
plot(x,e);
xlabel('x (m)');
ylabel('\epsilon');
title('End Profile (t = 10s)');

sgtitle('Forced Solution Example (Serre) \Delta x = 0.0469m')


