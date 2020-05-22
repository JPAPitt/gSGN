% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/ImpDisp/SmoothDB/alpha0p1/timeseries/';

figure;
begi = 00;
bege = 5;
for k = begi:bege
    if k < 10
        expwdir = strcat(wdir,'0',int2str(k), '/');
    else 
        expwdir = strcat(wdir,int2str(k), '/');
    end 

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

    x1 = xhGuend(:,1);
    h1 = xhGuend(:,2);
    G1 = xhGuend(:,3);
    u1 = xhGuend(:,4);

    x1a = xhGuendana(:,1);
    h1a = xhGuendana(:,2);
    G1a = xhGuendana(:,3);
    u1a = xhGuendana(:,4);
    
    plot(x1,h1,'-','DisplayName',strcat("\beta_1 = \beta_2 = ",num2str(beta1val,2)));
    
    if k == begi
       hold on
    end
end
legend('show');
hold off





