% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = './Validation/Recon/Test1/';


Norms = importdata(strcat(wdir, 'L2err.dat'));
dx = Norms(:,1);
Normh = Norms(:,2);
NormG = Norms(:,3);
Normu = Norms(:,4);

figure;
loglog(dx,Normh,'s b',dx,NormG,'o r',dx,Normu,'^ k' , dx, dx.^2,'-k',dx, 0.05*dx.^4,'-k', 'MarkerSize',8);
legend('h','G', 'u');
xlabel('\Delta x');
ylabel('L_2');
% 
% clc;
% close all;
% 
% Energys = importdata(strcat(wdir, 'Energy.dat'));
% dxE = Energys(:,1);
% hE = Energys(:,2);
% GE = Energys(:,3);
% uhE = Energys(:,4);
% HE = Energys(:,5);
% 
% figure;
% loglog(dxE,hE,'s b',dxE,GE,'o r',dxE, uhE,'^ k', dxE,HE,'+ b', dx, dx.^2,'-k', 'MarkerSize',8);
% legend('h','G', 'uh', '\mathcal{H}', 'Slope 2','Location','northwest');
% xlabel('\Delta x');
% ylabel('C');
% axis([10^-3 1 10^-16 10]);
% xticks([10^-3,10^-2,10^-1,10^0]);
% yticks([10^-16,10^-12,10^-8,10^-4,10^0]);



