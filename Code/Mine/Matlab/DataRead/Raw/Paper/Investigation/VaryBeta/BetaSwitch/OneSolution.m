% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/BetaVary/SmoothDBInvestigation/"
wdirend = 'OneOffTests/2to1/dbalpha0p1/Serre2SWWE/'
wdir = strcat(wdirbase,wdirend)

EndA = importdata(strcat(wdir,'EndAna.dat' ));
End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));

t = End(:,1);
x = End(:,2);
h = End(:,3);
G = End(:,4);
u = End(:,5);

tA = EndA(:,1);
xA = EndA(:,2);
hA = EndA(:,3);
GA = EndA(:,4);
uA = EndA(:,5);


figure;
plot(x,h,'-b',x,hA,'--k');
xlabel('x (m)');
ylabel('h (m)');

figure;
plot(x,u,'-r',x,uA,'--k');
xlabel('x (m)');
ylabel('u (m/s)');

figure;
plot(x,G,'-g',x,GA,'--k');
xlabel('x (m)');
ylabel('G (m/s)');

