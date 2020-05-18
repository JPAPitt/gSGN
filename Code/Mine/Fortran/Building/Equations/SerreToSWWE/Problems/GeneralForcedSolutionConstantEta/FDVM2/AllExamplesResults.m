% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir1 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta1/';
wdir0p5 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta0p5/';
wdir0 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta0/';

exnum = '06/';

xhGueps1 = importdata(strcat(wdir1,exnum , 'AllVarsEnd.dat' ));
xhGueps0p5 = importdata(strcat(wdir0p5,exnum , 'AllVarsEnd.dat' ));
xhGueps0 = importdata(strcat(wdir0,exnum , 'AllVarsEnd.dat' ));

xhGueps1Init = importdata(strcat(wdir1,exnum , 'InitVars.dat' ));
xhGueps0p5Init = importdata(strcat(wdir0p5,exnum , 'InitVars.dat' ));
xhGueps0Init = importdata(strcat(wdir0,exnum, 'InitVars.dat' ));

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

figure;
subplot(3, 4, 1);
plot(x0,h0);
xlabel('x (m)');
ylabel('height (m)');
title('t = 0s');

subplot(3, 4, 2);
plot(x,h);
xlabel('x (m)');
ylabel('height (m)');
title('Numerical Solution (t = 10s)');

subplot(3, 4, 3);
plot(x,hA);
xlabel('x (m)');
ylabel('height (m)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 4);
plot(x,((hA - h) ./ hA));
xlabel('x (m)');
ylabel(' ^{(h_{forced} - h_{computed})}/_{h_{forced}}');
title('Error (t = 10s)');

subplot(3, 4, 5);
plot(x0,u0);
xlabel('x (m)');
ylabel('u (m/s)');
title('t = 0s');

subplot(3, 4, 6);
plot(x,u);
xlabel('x (m)');
ylabel('u (m/s)');
title('Numerical Solution (t = 10s)');

subplot(3, 4,7);
plot(x,uA);
xlabel('x (m)');
ylabel('u (m/s)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 8);
plot(x,(uA - u));
xlabel('x (m)');
ylabel(' ^{(u_{forced} - u_{computed})}');
title('Error (t = 10s)');

subplot(3, 4, 9);
plot(x0,G0);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('t = 0s');

subplot(3, 4, 10);
plot(x,G);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Numerical Solution (t = 10s)');

subplot(3, 4,11);
plot(x,GA);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 12);
plot(x,(GA - G));
xlabel('x (m)');
ylabel(' ^{(G_{forced} - G_{computed})}');
title('Error (t = 10s)');

sgtitle('Forced Solution Example \epsilon = 1 (Serre) \Delta x = 0.0469m')


init = xhGueps0p5Init;
xhGu  = xhGueps0p5;

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

figure;
subplot(3, 4, 1);
plot(x0,h0);
xlabel('x (m)');
ylabel('height (m)');
title('t = 0s');

subplot(3, 4, 2);
plot(x,h);
xlabel('x (m)');
ylabel('height (m)');
title('Numerical Solution (t = 10s)');

subplot(3, 4, 3);
plot(x,hA);
xlabel('x (m)');
ylabel('height (m)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 4);
plot(x,((hA - h) ./ hA));
xlabel('x (m)');
ylabel(' ^{(h_{forced} - h_{computed})}/_{h_{forced}}');
title('Error (t = 10s)');

subplot(3, 4, 5);
plot(x0,u0);
xlabel('x (m)');
ylabel('u (m/s)');
title('t = 0s');

subplot(3, 4, 6);
plot(x,u);
xlabel('x (m)');
ylabel('u (m/s)');
title('Numerical Solution (t = 10s)');

subplot(3, 4,7);
plot(x,uA);
xlabel('x (m)');
ylabel('u (m/s)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 8);
plot(x,(uA - u));
xlabel('x (m)');
ylabel(' ^{(u_{forced} - u_{computed})}');
title('Error (t = 10s)');

subplot(3, 4, 9);
plot(x0,G0);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('t = 0s');

subplot(3, 4, 10);
plot(x,G);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Numerical Solution (t = 10s)');

subplot(3, 4,11);
plot(x,GA);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 12);
plot(x,(GA - G));
xlabel('x (m)');
ylabel(' ^{(G_{forced} - G_{computed})}');
title('Error (t = 10s)');

sgtitle('Forced Solution Example \epsilon = 0.5 \Delta x = 0.0469m')

init = xhGueps0Init;
xhGu  = xhGueps0;

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

figure;
subplot(3, 4, 1);
plot(x0,h0);
xlabel('x (m)');
ylabel('height (m)');
title('t = 0s');

subplot(3, 4, 2);
plot(x,h);
xlabel('x (m)');
ylabel('height (m)');
title('Numerical Solution (t = 10s)');

subplot(3, 4, 3);
plot(x,hA);
xlabel('x (m)');
ylabel('height (m)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 4);
plot(x,((hA - h) ./ hA));
xlabel('x (m)');
ylabel(' ^{(h_{forced} - h_{computed})}/_{h_{forced}}');
title('Error (t = 10s)');

subplot(3, 4, 5);
plot(x0,u0);
xlabel('x (m)');
ylabel('u (m/s)');
title('t = 0s');

subplot(3, 4, 6);
plot(x,u);
xlabel('x (m)');
ylabel('u (m/s)');
title('Numerical Solution (t = 10s)');

subplot(3, 4,7);
plot(x,uA);
xlabel('x (m)');
ylabel('u (m/s)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 8);
plot(x,(uA - u));
xlabel('x (m)');
ylabel(' ^{(u_{forced} - u_{computed})}');
title('Error (t = 10s)');

subplot(3, 4, 9);
plot(x0,G0);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('t = 0s');

subplot(3, 4, 10);
plot(x,G);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Numerical Solution (t = 10s)');

subplot(3, 4,11);
plot(x,GA);
xlabel('x (m)');
ylabel('G (m^2/s)');
title('Forced Solution (t = 10s)');

subplot(3, 4, 12);
plot(x,(GA - G));
xlabel('x (m)');
ylabel(' ^{(G_{forced} - G_{computed})}');
title('Error (t = 10s)');

sgtitle('Forced Solution Example \epsilon = 0 (SWWE) \Delta x = 0.0469m')




