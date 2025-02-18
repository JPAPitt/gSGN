% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data


wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/Forced/NV-beta1beta2Vary/';


OutEps = importdata(strcat(wdir, 'Norms.dat'));
dxseps = OutEps(:,1);
Normhseps = OutEps(:,2);
NormGseps = OutEps(:,3);
Normuseps = OutEps(:,4);

figure;
loglog(dxseps,Normhseps,'s b',dxseps,NormGseps,'o r',dxseps,Normuseps,'^ k','MarkerSize',8)
grid off
legend('h','G', 'u','Location','northwest')
xlabel('\Delta x')
ylabel('L_2')
axis([10^-4 10 10^-12 1]);
xticks([10^-4,10^-3,10^-2,10^-1,10^0,10]);
yticks([10^-12,10^-9,10^-6,10^-3,1]);
matlab2tikz('NormResults.tex');




