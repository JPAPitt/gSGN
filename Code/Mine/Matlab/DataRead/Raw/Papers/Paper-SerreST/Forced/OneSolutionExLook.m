% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SerreST/ForcedSolutionShort/07/"


End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));


t0 = Init(1,1);
x0 = Init(:,2);
h0 = Init(:,3);
G0 = Init(:,4);
u0 = Init(:,5);

t = End(1,1);
x = End(:,2);
h = End(:,3);
hA = End(:,4);
G = End(:,5);
GA = End(:,6);
u = End(:,7);
uA = End(:,8);


figure;
plot(x,h,'-r');
hold on;
plot(x,hA,'-b');
plot(x,h0,'-k');
hold off
xlabel('x (m)');
ylabel('h (m)');
legend('hide')



