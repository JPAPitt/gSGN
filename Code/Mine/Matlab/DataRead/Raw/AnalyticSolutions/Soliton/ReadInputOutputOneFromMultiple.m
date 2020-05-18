% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/GenSWWE/Forced/eps1o3/';
wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/LoopSoliton/';
exp = '06';
ptsep = 15;
linesep = 5;

initfile = strcat(wdir,'/',exp ,'/','InitVal.dat' );
endfile = strcat(wdir,'/',exp ,'/', 'EndVals.dat' );
endanafile = strcat(wdir,'/',exp ,'/', 'EndAnaVals.dat' );

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


figure;
plot(x0(1:linesep:end),h0(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('h (m)');
axis([-50 100 0.9 1.8]);
xticks([-50,0,50,100]);
yticks([1,1.2,1.4,1.6,1.8]);
legend('hide')
matlab2tikz('ExampleResultsHInit.tex');

clc;
close all;

figure;
plot(x1(1:linesep:end),h1a(1:linesep:end),'-b',x1(1:ptsep:end),h1(1:ptsep:end),'. k');
xlabel('x (m)');
ylabel('h (m)');
axis([-50 100 0.9 1.8]);
xticks([-50,0,50,100]);
yticks([1,1.2,1.4,1.6,1.8]);
legend('hide')
matlab2tikz('ExampleResultsHEnd.tex');


clc;
close all;

figure;
plot(x0(1:linesep:end),u0(1:linesep:end),'-r');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-50 100 -0.1 2]);
xticks([-50,0,50,100]);
yticks([0,0.5,1.0,1.5,2]);
legend('hide')
matlab2tikz('ExampleResultsUInit.tex');

clc;
close all;

figure;
plot(x1(1:linesep:end),u1a(1:linesep:end),'-r',x1(1:ptsep:end),u1(1:ptsep:end),'. k');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-50 100 -0.1 2]);
xticks([-50,0,50,100]);
yticks([0,0.5,1.0,1.5,2]);
legend('hide')
matlab2tikz('ExampleResultsUEnd.tex');

clc;
close all;
figure;
plot(x0(1:linesep:end),G0(1:linesep:end),'-g');
xlabel('x (m)');
ylabel('G (m^2/s)');
title('t = 0s');
axis([-50 100 -0.1 4]);
xticks([-50,0,50,100]);
yticks([0,1,2,3,4]);
legend('hide')
matlab2tikz('ExampleResultsGInit.tex');

clc;
close all;
figure;
plot(x1(1:linesep:end),G1a(1:linesep:end),'-g',x1(1:ptsep:end),G1(1:ptsep:end),'. k');
xlabel('x (m)');
ylabel('G (m^2/s)');
title('t = 10s');
axis([-50 100 -0.1 4]);
xticks([-50,0,50,100]);
yticks([0,1,2,3,4]);
legend('hide')
matlab2tikz('ExampleResultsGEnd.tex');
clc;
close all;


