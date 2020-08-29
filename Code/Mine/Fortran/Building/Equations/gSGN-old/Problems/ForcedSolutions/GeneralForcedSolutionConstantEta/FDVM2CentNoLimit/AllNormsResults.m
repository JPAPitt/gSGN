% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir1 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta1/';
wdir0p5 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta0p5/';
wdir0 = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/Serre2SWWE/Forced/eta0/';


OutEps1 = importdata(strcat(wdir1, 'Norms.dat'));
dxseps1 = OutEps1(:,1);
Normhseps1 = OutEps1(:,2);
NormGseps1 = OutEps1(:,3);
Normuseps1 = OutEps1(:,4);

figure;
slope2eps1 = 0.5*dxseps1.^(2);
loglog(dxseps1,Normhseps1,'s b',dxseps1,NormGseps1,'o r',dxseps1,Normuseps1,'^ k',dxseps1,slope2eps1 ,'-- k', 'MarkerSize',8)
grid on
legend('h','G', 'u', 'Slope 2','Location','northwest')
xlabel('\Delta x')
ylabel('L_2 Errors')
title('Forced Solution - Gaussian Epsilon = 1 (Serre)')
axis([10^-3 1 10^-8 10]);

OutEps0p5 = importdata(strcat(wdir0p5, 'Norms.dat'));
dxseps0p5 = OutEps0p5(:,1);
Normhseps0p5 = OutEps0p5(:,2);
NormGseps0p5 = OutEps0p5(:,3);
Normuseps0p5 = OutEps0p5(:,4);

figure;
slope2eps0p5 = 0.5*dxseps0p5.^(2);
loglog(dxseps0p5,Normhseps0p5,'s b',dxseps0p5,NormGseps0p5,'o r',dxseps0p5,Normuseps0p5,'^ k',dxseps0p5,slope2eps0p5 ,'-- k', 'MarkerSize',8)
grid on
legend('h','G', 'u', 'Slope 2','Location','northwest')
xlabel('\Delta x')
ylabel('L_2 Errors')
title('Forced Solution - Gaussian Epsilon = 0.5')
axis([10^-3 1 10^-8 10]);

OutEps0 = importdata(strcat(wdir0, 'Norms.dat'));
dxseps0 = OutEps0(:,1);
Normhseps0 = OutEps0(:,2);
NormGseps0 = OutEps0(:,3);
Normuseps0 = OutEps0(:,4);

figure;
slope2eps0 = 0.5*dxseps0.^(2);
loglog(dxseps0,Normhseps0,'s b',dxseps0,NormGseps0,'o r',dxseps0,Normuseps0,'^ k',dxseps0,slope2eps0 ,'-- k', 'MarkerSize',8)
grid on
legend('h','G', 'u', 'Slope 2','Location','northwest')
xlabel('\Delta x')
ylabel('L_2 Errors')
title('Forced Solution - Gaussian Epsilon = 0 (SWWE)')
axis([10^-3 1 10^-8 10]);




