% Process Fortran Outputs

% clc;
% clear all;
% close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';

%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaFunc/SWWE2Serre/SmoothDB/alpha0p1/timeseries/exp4/";
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/BetaFunc/Serre2SWWE/SmoothDB/alpha0p1/timeseries/exp6/";


%loop over time 
maxfile = 2000;
for k = 500:maxfile
    file = strcat(wdir,compose("%10d",k) ,'.dat' );

    xhGu = importdata(file);
    ctval = xhGu(1,1);
    betaval = xhGu(1,2);
    x1 = xhGu(:,3);
    h1 = xhGu(:,4);
    plot(x1,h1,'-b');
    title('Animated Plot of SWWE2Serre')
    xlabel('x');
    ylabel('h');
    legend(strcat('t=',num2str(ctval),'s' , '   \beta_1 = ',num2str(betaval)))
    drawnow;   
end 





