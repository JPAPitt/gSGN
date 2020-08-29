% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SerreST/ForcedSolution/01/"

dx = 3.1254883575558681E-002;

scenbeg = -100;
scenend = 100;

linestart = -100;
lineend = 100;
linestarti = 1+ round(((linestart - scenbeg)/dx));
lineendi = round(((lineend - scenbeg)/dx));
linesep = 1;
dotsep = 10;

End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));


t0 = Init(1,1);
x0 = Init(linestarti:linesep:lineendi,2);
h0 = Init(linestarti:linesep:lineendi,3);
G0 = Init(linestarti:linesep:lineendi,4);
u0 = Init(linestarti:linesep:lineendi,5);

t = End(1,1);
xdot = End(linestarti:dotsep:lineendi,2);
xline = End(linestarti:linesep:lineendi,2);
h = End(linestarti:dotsep:lineendi,3);
hA = End(linestarti:linesep:lineendi,4);
G = End(linestarti:dotsep:lineendi,5);
GA = End(linestarti:linesep:lineendi,6);
u = End(linestarti:dotsep:lineendi,7);
uA = End(linestarti:linesep:lineendi,8);


figure;
plot(xdot,h,'.r');
hold on;
plot(xline,hA,'-b');
plot(x0,h0,'-k');
hold off
xlabel('x (m)');
ylabel('h (m)');
axis([-25 75 0.9,1.6]);
xticks([-25,0,25,50,75]);
yticks([0.9,1,1.1,1.2,1.3,1.4,1.5,1.6]);
legend('hide')



figure;
plot(xdot,G,'.r');
hold on;
plot(xline,GA,'-b');
plot(x0,G0,'-k');
hold off
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-25 75 -0.1 0.5]);
xticks([-25,0,25,50,75]);
yticks([-0.1,0,0.1,0.2,0.3,0.4,0.5]);
legend('hide')


figure;
plot(xdot,u,'.r');
hold on;
plot(xline,uA,'-b');
plot(x0,u0,'-k');
hold off
xlabel('x (m)');
ylabel('u (m/s)');
axis([-25 75 -0.1 0.4]);
xticks([-25,0,25,50,75]);
yticks([-0.1,0,0.1,0.2,0.3,0.4]);
legend('hide');


