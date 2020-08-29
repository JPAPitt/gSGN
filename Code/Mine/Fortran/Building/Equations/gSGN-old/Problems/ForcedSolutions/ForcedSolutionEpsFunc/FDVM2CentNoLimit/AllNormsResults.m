% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/GenSWWE/ForcedEpsFunc/expTimeDeriv1/';


OutEps = importdata(strcat(wdir, 'Norms.dat'));
dxseps = OutEps(:,1);
Normhseps = OutEps(:,2);
NormGseps = OutEps(:,3);
Normuseps = OutEps(:,4);

figure;
slope2eps = 0.5*dxseps.^(2);
loglog(dxseps,Normhseps,'s b',dxseps,NormGseps,'o r',dxseps,Normuseps,'^ k',dxseps,slope2eps ,'-- k', 'MarkerSize',8)
grid on
legend('h','G', 'u', 'Slope 2','Location','northwest')
xlabel('\Delta x')
ylabel('L_2 Errors')
title('Forced Solution - Gaussian (Serre)')
axis([10^-3 1 10^-8 10]);




