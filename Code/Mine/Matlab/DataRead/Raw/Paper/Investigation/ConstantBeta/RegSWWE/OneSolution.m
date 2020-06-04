% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/AspRat2to1/DBalph0p1/"
subdir = 'Beta00/dx05/'

expwdir = strcat(wdir,subdir)
EndA = importdata(strcat(expwdir,'EndAna.dat' ));
End = importdata(strcat(expwdir,'End.dat' ));
% Init = importdata(strcat(wdir,'InitVal.dat' ));

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
subplot(1,2,1);
plot(x,h,'-b',x,hA,'--k');
xlabel('x (m)');
ylabel('h (m)');

subplot(1,2,2);
plot(x,h - hA,'-r');
xlabel('x (m)');
ylabel('h - hA (m)');

figure;
subplot(1,2,1);
plot(x,u,'-b',x,uA,'--k');
xlabel('x (m)');
ylabel('u (m/s)');

subplot(1,2,2);
plot(x,u - uA,'-r');
xlabel('x (m)');
ylabel('u - uA (m/s)');



figure;
subplot(1,2,1);
plot(x,G,'-b',x,GA,'--k');
xlabel('x (m)');
ylabel('G (m/s)');

subplot(1,2,2);
plot(x,G - GA,'-r');
xlabel('x (m)');
ylabel('G - GA (m)');
