% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/Serre2SWWECloser/";
wdir = "AspRat2to1/DBalpha0p1/Beta51";
dxdir = '/dx06/';
dx = 7.8126220722198776E-003;
scenbeg = -250;
scenend = 250;

u2 = 1.30584;
h2 = 1.45384;
ga = 9.81;


linestart = -200;
lineend = 200;
linestarti = round(((linestart - scenbeg)/dx));
lineendi = round(((lineend - scenbeg)/dx));
linesep = 10;

expdir = strcat(wdirbase,wdir,dxdir);

End = importdata(strcat(expdir,'End.dat' ));
param = fileread(strcat(expdir ,'Params.dat' ));

beta1str = extractBetween(param,"beta1 :","beta2");
beta1 = str2double(beta1str{1,1});

beta2str = extractAfter(param,"beta2 :");
beta2 = str2double(beta2str);

t = End(1,1);
x = End(linestarti:linesep:lineendi,2);
h = End(linestarti:linesep:lineendi,3);
G = End(linestarti:linesep:lineendi,4);
u = End(linestarti:linesep:lineendi,5);

Region1lb = -300;
Region1ub = (u2 - sqrt(ga*h2))*t;
Region2ub = (u2 - sqrt(ga*h2)*sqrt(beta2/ (2.0/3.0 + beta1)))*t;
Region3ub = (u2)*t;
Region4ub = (u2 + sqrt(ga*h2)*sqrt(beta2/ (2.0/3.0 + beta1)))*t;
Region5ub = (u2 + sqrt(ga*h2))*t;
Region6ub = 300;

ymin = 0.7;
ymax = 2.3;
x2 = [[Region1lb Region1lb],fliplr([Region1ub Region1ub])];
inBetween = [[ymin ymax], fliplr([ymin ymax])];
fill(x2, inBetween, 'b','facealpha',.1,'LineStyle','--');
hold on

x2 = [[Region1ub Region1ub],fliplr([Region2ub Region2ub])];
inBetween = [[ymin ymax], fliplr([ymin ymax])];
fill(x2, inBetween, 'r','facealpha',.1,'LineStyle','none');

x2 = [[Region2ub Region2ub],fliplr([Region3ub Region3ub])];
inBetween = [[ymin ymax], fliplr([ymin ymax])];
fill(x2, inBetween, 'g','facealpha',.1,'LineStyle','none');

x2 = [[Region3ub Region3ub],fliplr([Region4ub Region4ub])];
inBetween = [[ymin ymax], fliplr([ymin ymax])];
fill(x2, inBetween, 'g','facealpha',.1,'LineStyle','none');

x2 = [[Region4ub Region4ub],fliplr([Region5ub Region5ub])];
inBetween = [[ymin ymax], fliplr([ymin ymax])];
fill(x2, inBetween, 'r','facealpha',.1,'LineStyle','--');

x2 = [[Region5ub Region5ub],fliplr([Region6ub Region6ub])];
inBetween = [[ymin ymax], fliplr([ymin ymax])];
fill(x2, inBetween, 'b','facealpha',.1,'LineStyle','none');


plot(x,h,'-','DisplayName',strcat('\beta_2=',num2str(beta2)));
xlabel('x (m)');
ylabel('h (m)');
axis([-200 200 0.8 2.2]);
xticks([-200,-100,0,100,200]);
yticks([0.8,1,1.2,1.4,1.6,1.8,2,2.2]);
hold off
legend('hide')

matlab2tikz('hRegionsSerre.tex');
