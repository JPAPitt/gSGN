% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/ConstantBeta/SWWE/DB/timeseries/'
initfile = strcat(wdir, 'InitVal.dat' );
endfile = strcat(wdir, 'EndVals.dat' );
endanafile = strcat(wdir, 'EndAnaVals.dat' );
Energyfile = strcat(wdir, 'Energy.dat' );
Paramfile = strcat(wdir, 'Params.dat' );

ptsep = 100;
linesep = 25;


EnergyInfo = fileread(Energyfile);
InitEnergies = extractBetween(EnergyInfo,"H","End");
InitEnergies = str2num(InitEnergies{1,1});
EndEnergies = extractAfter(extractBetween(EnergyInfo,"End","Rel"),'H');
EndEnergies = str2num(EndEnergies{1,1});

ParamInfo = fileread(Paramfile);
EndTime = extractBetween(ParamInfo,"time :","dt");
EndTime = str2num(EndTime{1,1});

g = 9.81;
hl = 2;
hr = 1;
ChangeInTotal = [0,-g/2*(hr^2 - hl^2)*(EndTime),-g/2*(hr^2 - hl^2)*(EndTime),0];

writematrix(InitEnergies,'Energy.txt');
writematrix(EndEnergies,'Energy.txt','WriteMode','append');
writematrix(ChangeInTotal,'Energy.txt','WriteMode','append');

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

%Initial Conditions
% figure;
% plot(x0(1:linesep:end),h0(1:linesep:end),'-b');
% xlabel('x (m)');
% ylabel('h (m)');
% axis([-100 100 0.9 2.1]);
% xticks([-100,-50,0,50,100]);
% yticks([1,1.25,1.5,1.75,2]);
% legend('hide')
% matlab2tikz('DBInith.tex');
% 
% clc;
% close all;
% 
% figure;
% plot(x0(1:linesep:end),u0(1:linesep:end),'-r');
% xlabel('x (m)');
% ylabel('u (m/s)');
% axis([-100 100 -0.1 0.1]);
% xticks([-100,-50,0,50,100]);
% yticks([-0.1,0,0.1]);
% legend('hide')
% matlab2tikz('DBInitu.tex');
% 
% clc;
% close all;
% 
% figure;
% plot(x0(1:linesep:end),G0(1:linesep:end),'-g');
% xlabel('x (m)');
% ylabel('G (m^2/s)');
% axis([-100 100 -0.1 0.1]);
% xticks([-100,-50,0,50,100]);
% yticks([-0.1,0,0.1]);
% legend('hide')
% matlab2tikz('DBInitG.tex');

%Final Solution

figure;
plot(x1(1:linesep:end),h1a(1:linesep:end),'-b',x1(1:ptsep:end),h1(1:ptsep:end),'.k');
xlabel('x (m)');
ylabel('h (m)');
axis([-100 100 0.9 2.1]);
xticks([-100,-50,0,50,100]);
yticks([1,1.25,1.5,1.75,2]);
legend('hide')
matlab2tikz('DBEndh.tex');
 
clc;
close all;

figure;
plot(x1(1:linesep:end),u1a(1:linesep:end),'-r',x1(1:ptsep:end),u1(1:ptsep:end),'.k');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-100 100 -0.1 1.5]);
xticks([-100,-50,0,50,100]);
yticks([0,0.5,1.0,1.5]);
legend('hide')
matlab2tikz('DBEndu.tex');
 
clc;
close all;

figure;
plot(x1(1:linesep:end),G1a(1:linesep:end),'-g',x1(1:ptsep:end),G1(1:ptsep:end),'.k');
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-100 100 -0.1 2]);
xticks([-100,-50,0,50,100]);
yticks([0,0.5,1.0,1.5,2]);
legend('hide')
matlab2tikz('DBEndG.tex');
 
clc;
close all;




