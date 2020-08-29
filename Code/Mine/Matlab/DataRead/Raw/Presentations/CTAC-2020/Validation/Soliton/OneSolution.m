% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/SerreSoliton/06/"


EndA = importdata(strcat(wdir,'EndAna.dat' ));
End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));

skip = 2;
t = End(1,1);
x = End(1:skip:end,2);
h = End(1:skip:end,3);
G = End(1:skip:end,4);
u = End(1:skip:end,5);

tA = EndA(1,1);
xA = EndA(:,2);
hA = EndA(:,3);
GA = EndA(:,4);
uA = EndA(:,5);

ti = Init(1,1);
xi = Init(:,2);
hi = Init(:,3);
Gi = Init(:,4);
ui = Init(:,5);

figure;
plot(xi,hi,'-k',xA,hA,'-b');
xlabel('x (m)');
ylabel('h (m)');
xlim([-50 150]);
ylim([1 1.8]);
xticks([-50,0,50,100,150]);
yticks([1,1.2,1.4,1.6,1.8]);
legend('t=0','Analytic t = 30s')
matlab2tikz('Solution.tex');
close all;

figure;
plot(xi,hi,'-k',xA,hA,'-b',x,h,'.r');
xlabel('x (m)');
ylabel('h (m)');
xlim([-50 150]);
ylim([1 1.8]);
xticks([-50,0,50,100,150]);
yticks([1,1.2,1.4,1.6,1.8]);
legend('t=0','Analytic t = 30s', 'Numeric t = 30s')
matlab2tikz('Numeric.tex');
close all;


