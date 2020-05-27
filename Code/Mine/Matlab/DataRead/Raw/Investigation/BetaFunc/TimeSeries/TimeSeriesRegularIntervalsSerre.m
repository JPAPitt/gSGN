% Process Fortran Outputs

% clc;
% clear all;
% close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaConstant/Serre/SmoothDB/1to0p1/alpha0p1/timeseries/exp6/"
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/Beta1LimitOpt1/SmoothDB/alpha0p1/timeseries/exp6/"

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/Beta1LimitOpt1/SmoothDB/2to1/alpha0p1/timeseries/exp1/"

linesep = 1;

file = strcat(wdir, 'InitVal.dat' );
filedat = importdata(file);
x = filedat(:,1);
h = filedat(:,2);
G = filedat(:,3);
u = filedat(:,4);

figure;
plot(x(1:linesep:end),h(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('h (m)');
title('t=0s');
legend('hide');


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
%times = [0.206,0.207,0.2071,0.2072,0.2073,0.2074,0.2075,0.2077,0.208]
times = [1,10,15]
for k = 1:size(times,2)
    currenttime = times(k);
    currenttimefileindex = round(currenttime / filedt);
    
    filei = strcat(wdir,compose("%10d",currenttimefileindex) ,'.dat' );
    fileidat = importdata(filei);
    fileitime = fileidat(1,1);
    fileip1 = strcat(wdir,compose("%10d",currenttimefileindex + 1) ,'.dat' );
    fileip1dat = importdata(fileip1);
    fileip1time = fileip1dat(1,1);
    
    filedat = fileidat;
    if (abs(fileip1time - currenttime) < abs(fileitime - currenttime) )
        filedat = fileip1dat;
    end
    
    x = filedat(:,2);
    h = filedat(:,4);
%     G = filedat(:,4);
%     u = filedat(:,5);
    
    figure;
    plot(x,h,'-b');
    title(strcat('t=',num2str(currenttime),'s' , '   \beta_1 = 0'));
    legend('hide');
    
    
end 



file = strcat(wdir, 'EndVals.dat' );
filedat = importdata(file);
x = filedat(:,1);
h = filedat(:,2);
% G = filedat(:,3);
% u = filedat(:,4);

figure;
plot(x(1:linesep:end),h(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('h (m)');
title('t=60s');
legend('hide');






