% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SerreST/PeakonUp-10s/"


OutEps = importdata(strcat(wdir, 'Norms.dat'));
dxseps = OutEps(:,1);
Normhseps = OutEps(:,2);
NormGseps = OutEps(:,3);
Normuseps = OutEps(:,4);
H1hseps = OutEps(:,5);
H1Gseps = OutEps(:,6);
H1useps = OutEps(:,7);

figure;
loglog(dxseps,Normhseps,'s b',dxseps,NormGseps,'o r',dxseps,Normuseps,'^ k' , dxseps, dxseps.^2, '-', dxseps, dxseps, '-')
grid off
legend('h','G', 'u','Location','northwest')
xlabel('\Delta x')
ylabel('L_2')
legend('show')



figure;
loglog(dxseps,H1hseps,'s b',dxseps,H1Gseps,'o r',dxseps,H1useps,'^ k' , dxseps, dxseps.^2, '-', dxseps, dxseps, '-')
grid off
legend('h','G', 'u','Location','northwest')
xlabel('\Delta x')
ylabel('H_1')
legend('show')




