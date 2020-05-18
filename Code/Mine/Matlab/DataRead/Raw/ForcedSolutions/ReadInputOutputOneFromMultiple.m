% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/Forced/NV-beta1beta2Vary/';
exnum = '06/';
ptsep = 15;
linesep = 5;

xhGu = importdata(strcat(wdir,exnum , 'OutVars.dat' ));
init = importdata(strcat(wdir,exnum , 'InitVars.dat' ));

x = xhGu(:,1);
hA = xhGu(:,2);
h = xhGu(:,3);
GA = xhGu(:,4);
G = xhGu(:,5);
uA = xhGu(:,6);
u = xhGu(:,7);
beta1 = xhGu(:,8);
beta2 = xhGu(:,9);

x0 = init(:,1);
h0 = init(:,2);
G0 = init(:,3);
u0 = init(:,4);
beta1t0 = init(:,5);
beta2t0 = init(:,6);

figure;
plot(x0(1:linesep:end),h0(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('h (m)');
axis([-50 100 0.9 1.6]);
xticks([-50,0,50,100]);
yticks([1,1.2,1.4,1.6]);
legend('hide')
matlab2tikz('ExampleResultsHInit.tex');

clc;
close all;
figure;
plot(x(1:linesep:end),hA(1:linesep:end),'-b',x(1:ptsep:end),h(1:ptsep:end),'. k');
xlabel('x (m)');
ylabel('h (m)');
axis([-50 100 0.9 1.6]);
xticks([-50,0,50,100]);
yticks([1,1.2,1.4,1.6]);
legend('hide')
matlab2tikz('ExampleResultsHEnd.tex');

clc;
close all;
figure;
plot(x0(1:linesep:end),u0(1:linesep:end),'-r');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-50 100 -0.1 0.4]);
xticks([-50,0,50,100]);
yticks([-0.1,0,0.1,0.2,0.3,0.4]);
legend('hide')
matlab2tikz('ExampleResultsUInit.tex');

clc;
close all;
figure;
plot(x(1:linesep:end),uA(1:linesep:end),'-r',x(1:ptsep:end),u(1:ptsep:end),'. k');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-50 100 -0.1 0.4]);
xticks([-50,0,50,100]);
yticks([-0.1,0,0.1,0.2,0.3,0.4]);
legend('hide')
matlab2tikz('ExampleResultsUEnd.tex');

clc;
close all;
figure;
plot(x0(1:linesep:end),u0(1:linesep:end),'-g');
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-50 100 -0.1 0.4]);
xticks([-50,0,50,100]);
yticks([-0.1,0,0.1,0.2,0.3,0.4]);
legend('hide')
matlab2tikz('ExampleResultsGInit.tex');

clc;
close all;
figure;
plot(x(1:linesep:end),uA(1:linesep:end),'-g',x(1:ptsep:end),u(1:ptsep:end),'. k');
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-50 100 -0.1 0.4]);
xticks([-50,0,50,100]);
yticks([-0.1,0,0.1,0.2,0.3,0.4]);
legend('hide')
matlab2tikz('ExampleResultsGEnd.tex');

clc;
close all;
figure;
plot(x(1:linesep:end),beta1(1:linesep:end),'-b',x0(1:linesep:end),beta1t0(1:linesep:end),'-k');
xlabel('x (m)');
ylabel('\beta_1');
axis([-50 100 -1 2]);
xticks([-50,0,50,100]);
yticks([-1,-0.5,0,0.5,1,1.5,2]);
legend('hide')
matlab2tikz('ExampleResultsBeta1.tex');

clc;
close all;
figure;
plot(x(1:linesep:end),beta2(1:linesep:end),'-b',x0(1:linesep:end),beta2t0(1:linesep:end),'-k');
xlabel('x (m)');
ylabel('\beta_2');
axis([-50 100 -0.1 3]);
xticks([-50,0,50,100]);
yticks([0,1,2,3]);
legend('hide')
matlab2tikz('ExampleResultsBeta2.tex');


