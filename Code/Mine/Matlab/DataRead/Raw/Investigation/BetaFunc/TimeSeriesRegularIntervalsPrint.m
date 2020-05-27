% Process Fortran Outputs

% clc;
% clear all;
% close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaFunc/SWWE2Serre/SmoothDB/alpha0p1/timeseries/exp2/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaFunc/Serre2SWWE/SmoothDB/alpha0p1/timeseries/exp4/";

%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaConstant/Serre/SmoothDB/alpha0p1/timeseries/exp1/"

linesep = 20;

file = strcat(wdir, 'InitVal.dat' );
filedat = importdata(file);
x = filedat(:,1);
h = filedat(:,2);
G = filedat(:,3);
u = filedat(:,4);

figure;
plot(x,h,'-b');
xlabel('x (m)');
ylabel('h (m)');
axis([-300 300 0.9 2.1]);
xticks([-300,-150,0,150,300]);
yticks([1,1.25,1.5,1.75,2]);
title('t=0s');

matlab2tikz('SWWE2Serre10to11A.tex');

clc;
close all;

paramfile = strcat(wdir ,'Params.dat' );
param = fileread(paramfile);
    
dtstr = extractBetween(param,"dx*(dt/dx)  :","gravity");
dt = str2double(dtstr{1,1});

file0 = strcat(wdir,compose("%10d",0) ,'.dat' );
file0dat = importdata(file0);
startt = file0dat(1,1);
file1 = strcat(wdir,compose("%10d",1) ,'.dat' );
file1dat = importdata(file1);
nextt = file1dat(1,1);

filedt = nextt - startt;
times = [10,20,30,40,50]

for k = 1:size(times,2)
    currenttime = times(k);
    currenttimefileindex = round(currenttime / filedt);
    
    filei = strcat(wdir,compose("%10d",currenttimefileindex) ,'.dat' );
    fileidat = importdata(filei);
    fileitime = fileidat(1,1);
    fileip1 = strcat(wdir,compose("%10d",currenttimefileindex + 1) ,'.dat' );
    fileip1dat = importdata(fileip1);
    fileip1time = fileip1dat(1,1);
    
    
    if (abs(fileip1time - currenttime) < abs(fileitime - currenttime) )
        filedat = fileip1dat;
    else
        filedat = fileidat;
    end 
    
    beta1 = filedat(1,2);
    x = filedat(:,3);
    h = filedat(:,4);
    G = filedat(:,5);
    u = filedat(:,6);
    
    figure;
    plot(x(1:linesep:end),h(1:linesep:end),'-b');
    xlabel('x (m)');
    ylabel('h (m)');
    axis([-300 300 0.9 2.1]);
    xticks([-300,-150,0,150,300]);
    yticks([1,1.25,1.5,1.75,2]);
    title(strcat('t=',num2str(currenttime),'s' , '   \beta_1 = ',num2str(beta1)));
    
    matlab2tikz( strcat('SWWE2Serre10to11',65+ k,'.tex'));

    clc;
    close all;
    
end 



file = strcat(wdir, 'EndVals.dat' );
filedat = importdata(file);
x = filedat(:,1);
h = filedat(:,2);
G = filedat(:,3);
u = filedat(:,4);

figure;
plot(x,h,'-b');
xlabel('x (m)');
ylabel('h (m)');
axis([-300 300 0.9 2.1]);
xticks([-300,-150,0,150,300]);
yticks([1,1.25,1.5,1.75,2]);
title('t=60s');

matlab2tikz( strcat('SWWE2Serre10to11ht',65+ k + 1,'.tex'));

clc;
close all;




