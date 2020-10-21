% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "./Validation/Recon/SINE/"


OutEps = importdata(strcat(wdir, 'Norms.dat'));
dxseps = OutEps(:,1);
Normhmseps = OutEps(:,2);
Normhpseps = OutEps(:,3);

figure;
loglog(dxseps,Normhmseps,'s b',dxseps,Normhpseps,'^ r',dxseps,0.0001*dxseps.^2,'-',dxseps,0.000001*dxseps.^5,'-')
grid off
legend('hide')
xlabel('\Delta x')
ylabel('L_2')
% axis([10^-4 10 10^-8 10^2]);
% xticks([10^-4,10^-3,10^-2,10^-1,10^0,10]);
% yticks([10^-8,10^-6,10^-4,10^-2,1,10^2]);
% matlab2tikz('NormResults.tex');




