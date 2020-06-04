% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/BetaFunc/Beta1Beta2Switch/Serre2ImpSerre/SmoothDB/1to0p01/alpha0p1/timeseries/exp1/"

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

% matlab2tikz('Serre2SWWE20A.tex');
% 
% clc;
% close all;

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
times = [5,10,20,30];

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
    subplot(1,2,1)
    plot(x(1:linesep:end),h(1:linesep:end),'-b');
    title('h')
    subplot(1,2,2)
    plot(x(1:linesep:end),u(1:linesep:end),'-r');
    title('u')
    sgtitle(strcat('t=',num2str(currenttime),'s' , '   \beta_1 = ',num2str(beta1)));
    legend('hide');
    
%     matlab2tikz( strcat('Serre2SWWE20',65+ k,'.tex'));
% 
%     clc;
%     close all;
    
end 



file = strcat(wdir, 'EndVals.dat' );
filedat = importdata(file);
x = filedat(:,1);
h = filedat(:,2);
G = filedat(:,3);
u = filedat(:,4);

figure;
plot(x(1:linesep:end),h(1:linesep:end),'-b');
xlabel('x (m)');
ylabel('h (m)');
title('t=60s');
legend('hide');

% matlab2tikz( strcat('Serre2SWWE20',65+ k + 1,'.tex'));





