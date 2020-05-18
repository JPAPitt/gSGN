% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/SWWE/LoopDB/'


Norms = importdata(strcat(wdir, 'Norms.dat'));
dx = Norms(:,1);
Normh = Norms(:,2);
NormG = Norms(:,3);
Normu = Norms(:,4);

figure;
slope2 = dx.^(2);
loglog(dx,Normh,'s b',dx,NormG,'o r',dx,Normu,'^ k',dx,slope2 ,'-- k', 'MarkerSize',8)
grid on
legend('h','G', 'u', 'Slope 2','Location','northwest')
xlabel('\Delta x')
ylabel('L_2 Errors')
title('Soliton Solution (gSGN \beta_1=\beta_2 = 0 - Serre)')
axis([10^-3 1 10^-12 10]);


Energys = importdata(strcat(wdir, 'Energy.dat'));
dxE = Energys(:,1);
hE = Energys(:,2);
GE = Energys(:,3);
uhE = Energys(:,4);
HE = Energys(:,5);

figure;
slope2E = 0.5*dxE.^(2);
loglog(dxE,hE,'s b',dxE,GE,'o r',dxE, uhE,'^ k', dxE,HE,'+ b' ,dxE,slope2E ,'-- k', 'MarkerSize',8)
grid on
legend('h','G', 'uh', 'H', 'Slope 2','Location','northwest')
xlabel('\Delta x')
ylabel('Conservation Errors')
title('Soliton Solution (gSGN \beta_1=\beta_2 = 0 - Serre)')
axis([10^-3 1 10^-14 10]);



