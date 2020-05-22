% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/ImpDisp/Gaussian/timeseries/';

figure;
begi = 25;
bege = 30;
for k = begi:bege
    if k < 10
        expwdir = strcat(wdir,'0',int2str(k), '/');
    else 
        expwdir = strcat(wdir,int2str(k), '/');
    end 

    endfile = strcat(expwdir, 'EndVals.dat' );

    xhGuend = importdata(endfile);

    x1 = xhGuend(:,1);
    h1 = xhGuend(:,2);
    G1 = xhGuend(:,3);
    u1 = xhGuend(:,4);
    
    plot(x1,h1,'-','DisplayName',strcat("\beta = ",num2str(k/30.0) ));
    if k == begi
       hold on
    end
end
legend('show');
hold off





