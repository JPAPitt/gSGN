% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"


beta0 = importdata('REGSWWE_Beta00_ETS.dat');
beta0p5 = importdata('REGSWWE_Beta05_ETS.dat');
SWWE = importdata('SWWESolver_DB_ETS.dat');
SWWEDBS = importdata('SWWESolver_DBa0p1_ETS.dat');

tbeta0 = beta0(:,1);
Hbeta0 = beta0(:,5);

tbeta0p5 = beta0p5(:,1);
Hbeta0p5 = beta0p5(:,5);


plot(tbeta0,Hbeta0,'.','DisplayName','beta_2 = 0');
hold on;
plot(tbeta0p5,Hbeta0p5,'.','DisplayName','beta_2 = 0.5');
plot(SWWE(:,1),SWWE(:,5),'.','DisplayName','SWWE Solver DB');
plot(SWWEDBS(:,1),SWWEDBS(:,5),'.','DisplayName','SWWE Solver Smooth DB 0.1');

xlabel('t (s)');
ylabel('H (m^3/s^2)');

axis([0 35 6105 6145])
xticks([0,5,10,15,20,25,30,35])
yticks([6105,6110,6115,6120,6125,6130,6135,6140,6145])
legend('show')
hold off
% matlab2tikz('EnergyOverTime.tex');