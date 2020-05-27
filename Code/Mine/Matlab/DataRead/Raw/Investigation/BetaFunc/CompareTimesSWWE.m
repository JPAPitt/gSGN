% Process Fortran Outputs

% clc;
% clear all;
% close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaFunc/SWWE2Serre/SmoothDB/alpha0p1/timeseries/exp4/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaFunc/Serre2SWWE/SmoothDB/alpha0p1/timeseries/exp4/";

%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaConstant/Serre/SmoothDB/alpha0p1/timeseries/exp1/"
file = strcat(wdir, 'InitVal.dat' );

xhGuinit = importdata(file);
x0 = xhGuinit(:,1);
h0 = xhGuinit(:,2);
G0 = xhGuinit(:,3);
u0 = xhGuinit(:,4);

figure;
subplot(4, 4, 1);
plot(x0,h0,'-b');
title('t=0s');

%loop over time 
maxfile = 300;
sep = round(355/16);
for k = 2:15
    filenum =  k*sep;
    file = strcat(wdir,compose("%10d",filenum) ,'.dat' );
    
    xhGu = importdata(file);
    ctval = xhGu(1,1);
    betaval = xhGu(1,2);
    x1 = xhGu(:,3);
    h1 = xhGu(:,4);
    G1 = xhGu(:,5);
    u1 = xhGu(:,6);
    subplot(4, 4, k);
    plot(x1,h1,'-b');
    title(strcat('t=',num2str(ctval),'s' , '   \beta_1 = ',num2str(betaval)));

    
%     ctval = xhGu(1,1);
%     x1 = xhGu(:,2);
%     h1 = xhGu(:,3);
%     G1 = xhGu(:,4);
%     u1 = xhGu(:,5);
%     subplot(4, 4, k);
%     plot(x1,h1,'-b');
%     title(strcat('t=',num2str(ctval),'s'));

    
    
    
end 


endfile = strcat(wdir, 'EndVals.dat' );
endanafile = strcat(wdir, 'EndAnaVals.dat' );

xhGuend = importdata(endfile);
xhGuendana = importdata(endanafile);

x1 = xhGuend(:,1);
h1 = xhGuend(:,2);
G1 = xhGuend(:,3);
u1 = xhGuend(:,4);

x1a = xhGuendana(:,1);
h1a = xhGuendana(:,2);
G1a = xhGuendana(:,3);
u1a = xhGuendana(:,4);

subplot(4, 4, 16);
plot(x1,h1,'-b',x1,h1a,'--k');
title('t=60s');




