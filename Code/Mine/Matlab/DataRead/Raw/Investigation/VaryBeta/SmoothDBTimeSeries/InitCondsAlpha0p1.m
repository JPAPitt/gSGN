% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';
linesep = 25;
linesep1 = 1;
wdirSerre = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/ImpDisp/SmoothDB/alpha0p1/timeseries/00/';

initfile = strcat(wdirSerre, 'InitVal.dat' );


%make initial condition figure
xhGuinit = importdata(initfile);
x0 = xhGuinit(:,1);
h0 = xhGuinit(:,2);
G0 = xhGuinit(:,3);
u0 = xhGuinit(:,4);

figure;
plot(x0(1:linesep:end),h0(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('h (m)');
axis([-100 100 0.9 2.1]);
xticks([-100,-50,0,50,100]);
yticks([1,1.25,1.5,1.75,2]);
legend('hide')
matlab2tikz('SDBa0p1Inith.tex');

clc;
close all;

figure;
plot(x0(4500:linesep1:5500),h0(4500:linesep1:5500),'-b');
xlabel('x (m)');
ylabel('h (m)');
axis([-5 5 0.9 2.1]);
xticks([-5,-2.5,0,2.5,5]);
yticks([1,1.25,1.5,1.75,2]);
legend('hide')
matlab2tikz('SDBa0p1Inithz.tex');

clc;
close all;

figure;
plot(x0(1:linesep:end),u0(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-100 100 -0.1 0.1]);
xticks([-100,-50,0,50,100]);
yticks([-0.1, 0, 0.1]);
legend('hide')
matlab2tikz('SDBa0p1Initu.tex');

clc;
close all;

figure;
plot(x0(1:linesep:end),G0(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-100 100 -0.1 0.1]);
xticks([-100,-50,0,50,100]);
yticks([-0.1, 0, 0.1]);
legend('hide')
matlab2tikz('SDBa0p1InitG.tex');

clc;
close all;