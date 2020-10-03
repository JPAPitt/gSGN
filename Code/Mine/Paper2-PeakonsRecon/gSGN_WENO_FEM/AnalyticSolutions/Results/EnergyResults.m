% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Run/30s/"


OutEps = importdata(strcat(wdir, 'Energy.dat'));
dxseps = OutEps(:,1);
mass = OutEps(:,2);
G = OutEps(:,3);
mome = OutEps(:,4);
energy = OutEps(:,5);

figure;
loglog(dxseps,mass,'s b',dxseps,G,'o r',dxseps,mome,'^ k',dxseps,energy,'x g',dxseps,dxseps.^2,'-',dxseps,dxseps.^4,'-')
grid off
legend('hide')
xlabel('\Delta x')
ylabel('C_1')
% axis([10^-4 10 10^-8 10^2]);
% xticks([10^-4,10^-3,10^-2,10^-1,10^0,10]);
% yticks([10^-8,10^-6,10^-4,10^-2,1,10^2]);
% matlab2tikz('NormResults.tex');




