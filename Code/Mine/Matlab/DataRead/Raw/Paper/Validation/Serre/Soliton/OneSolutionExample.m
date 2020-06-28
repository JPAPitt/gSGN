% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/Validation/AnalyticSolutions/SerreSolitonWider/06/"

param = fileread(strcat(wdir ,'Params.dat' ));
dxstr = extractBetween(param,"dx","tstart");
dx = str2double(dxstr{1,1});

scenbeg = -200;
scenend = 200;

linestart = -25;
lineend = 150;
linestarti = round(((linestart - scenbeg)/dx));
lineendi = round(((lineend - scenbeg)/dx));
linesep = 1;
dotsep = 5;


EndA = importdata(strcat(wdir,'EndAna.dat' ));
End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));


t0 = Init(1,1);
x0 = Init(linestarti:linesep:lineendi,2);
h0 = Init(linestarti:linesep:lineendi,3);
G0 = Init(linestarti:linesep:lineendi,4);
u0 = Init(linestarti:linesep:lineendi,5);

t = End(1,1);
x = End(linestarti:dotsep:lineendi,2);
h = End(linestarti:dotsep:lineendi,3);
G = End(linestarti:dotsep:lineendi,4);
u = End(linestarti:dotsep:lineendi,5);

tA = EndA(1,1);
xA = EndA(linestarti:linesep:lineendi,2);
hA = EndA(linestarti:linesep:lineendi,3);
GA = EndA(linestarti:linesep:lineendi,4);
uA = EndA(linestarti:linesep:lineendi,5);


figure;
plot(x,h,'.r');
hold on;
plot(xA,hA,'-b');
plot(x0,h0,'-k');
hold off
xlabel('x (m)');
ylabel('h (m)');
axis([-25 150 0.9,1.8]);
xticks([-25,0,25,50,75,100,125,150]);
yticks([0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8]);
legend('hide')
matlab2tikz('h.tex');
close all;


figure;
plot(x,G,'.r');
hold on;
plot(xA,GA,'-b');
plot(x0,G0,'-k');
hold off
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-25 150 -0.5 4]);
xticks([-25,0,25,50,75,100,125,150]);
yticks([-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4]);
legend('hide')
matlab2tikz('G.tex');
close all;

figure;
plot(x,u,'.r');
hold on;
plot(xA,uA,'-b');
plot(x0,u0,'-k');
hold off
xlabel('x (m)');
ylabel('u (m/s)');
axis([-25 150 -0.2 1.8]);
xticks([-25,0,25,50,75,100,125,150]);
yticks([-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8]);
legend('hide');
matlab2tikz('u.tex');
close all;

