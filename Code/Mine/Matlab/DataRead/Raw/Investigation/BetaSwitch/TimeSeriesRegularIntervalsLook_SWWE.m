% Process Fortran Outputs

clc;
clear all;
close all;



% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/BetaFunc/BetaConstant/SWWE/SmoothDB/1to0p001/alpha0p5/timeseries/exp1/"

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
axis([-200 200 0.0 1.1]);

% matlab2tikz('SWWE0p1hA.tex');
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
    plot(x(1:linesep:end),h(1:linesep:end),'-b');
    title(strcat('h : ','t=',num2str(currenttime),'s' , '   \beta_1 = ',num2str(beta1)));
    axis([-200 200 0.0 1.1]);
    xlabel('x (m)');
    ylabel('h (m)');
   
    figure;
    plot(x(1:linesep:end),u(1:linesep:end),'-r');
    title(strcat('u : ','t=',num2str(currenttime),'s' , '   \beta_1 = ',num2str(beta1)));
    axis([-200 200 0.0 3]);
    xlabel('x (m)');
    ylabel('u (m/s)');
        
end 

file = strcat(wdir, 'EndVals.dat' );
filedat = importdata(file);
x = filedat(:,1);
h = filedat(:,2);
G = filedat(:,3);
u = filedat(:,4);

figure;
plot(x(1:linesep:end),h(1:linesep:end),'-b');
title(strcat('Numerical h : ','t=40s' , '   \beta_1 = ',num2str(beta1)));
axis([-200 200 0.0 1.1]);
xlabel('x (m)');
ylabel('h (m)');

figure;
plot(x(1:linesep:end),u(1:linesep:end),'-r');
title(strcat('Numerical u : ','t=40s' , '   \beta_1 = ',num2str(beta1)));
axis([-200 200 0.0 3]);
xlabel('x (m)');
ylabel('u (m/s)');


file = strcat(wdir, 'EndAnaVals.dat' );
filedat = importdata(file);
xA = filedat(:,1);
hA = filedat(:,2);
GA = filedat(:,3);
uA = filedat(:,4);

figure;
plot(x(1:linesep:end),h(1:linesep:end),'-b',xA(1:linesep:end),hA(1:linesep:end),'--k' );
title(strcat('Numerical and Analytical h : ','t=40s' , '   \beta_1 = ',num2str(beta1)));
axis([-200 200 0.0 1.1]);
xlabel('x (m)');
ylabel('h (m)');

figure;
plot(x(1:linesep:end),u(1:linesep:end),'-r',xA(1:linesep:end),uA(1:linesep:end),'--k' );
title(strcat('Numerical and Analytical h : ','t=40s' , '   \beta_1 = ',num2str(beta1)));
axis([-200 200 0.0 3]);
xlabel('x (m)');
ylabel('u (m/s)');






