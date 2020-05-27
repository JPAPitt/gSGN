% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

expwdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/ImpDisp/SmoothDB/alpha2/timeseries/00/';

figure;

linesep = 20;
u2 = 1.30584;
h2 = 1.45384;
ga = 9.81;


endfile = strcat(expwdir, 'EndVals.dat' );
endanafile = strcat(expwdir, 'EndAnaVals.dat' );
paramfile = strcat(expwdir, 'Params.dat' );

xhGuend = importdata(endfile);
xhGuendana = importdata(endanafile);
param = fileread(paramfile);

beta1 = extractBetween(param,"beta1 : ","beta2");
beta1val = str2double(beta1{1,1});
beta2 = extractAfter(param,"beta2 :");
beta2val = str2double(beta2);

endtime = extractBetween(param,"time :","dt");
endtimeval = str2double(endtime);

x1 = xhGuend(:,1);
h1 = xhGuend(:,2);
u1 = xhGuend(:,4);


Region1lb = -100;
Region1ub = (u2 - sqrt(ga*h2))*endtimeval;
Region2ub = (u2 - sqrt(ga*h2)*sqrt(beta2val/ (2.0/3.0 + beta1val)))*endtimeval;
Region3ub = (u2)*endtimeval;
Region4ub = (u2 + sqrt(ga*h2)*sqrt(beta2val/ (2.0/3.0 + beta1val)))*endtimeval;
Region5ub = (u2 + sqrt(ga*h2))*endtimeval;
Region6ub = 100;

ymin = 0.8
ymax = 2.2
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

plot(x1(1:linesep:end),h1(1:linesep:end),'-b','DisplayName',strcat("\beta_1 = ",num2str(beta1val,2), "  \beta_2 = ",num2str(beta2val,2)));

hold off;

legend('hide');
xlabel('x (m)');
ylabel('h (m)');
axis([-100 100 0.8 2.2]);
xticks([-100,-50,0,50,100]);
yticks([0.8,1,1.2,1.4,1.6,1.8,2,2.2]);

figstr = strcat('RegionhPlotBetaZero.tex')
matlab2tikz(figstr);

clc;
close all;

figure;
ymin = -1
ymax = 3
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

plot(x1(1:linesep:end),u1(1:linesep:end),'-b','DisplayName',strcat("\beta_1 = ",num2str(beta1val,2), "  \beta_2 = ",num2str(beta2val,2)));

hold off;

legend('hide');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-100 100 -1 3]);
xticks([-100,-50,0,50,100]);
yticks([-1,0,1,2,3]);

figstr = strcat('RegionuPlotBetaZero.tex')
matlab2tikz(figstr);










