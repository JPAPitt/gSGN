% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/RegSWWE/SmoothDB/alpha2/timeseries/';

figure;
begi = 0;
bege = 4;
for k = begi:bege
    if k < 10
        expwdir = strcat(wdir,'0',int2str(k), '/');
    else 
        expwdir = strcat(wdir,int2str(k), '/');
    end 

    endfile = strcat(expwdir, 'EndVals.dat' );
    endanafile = strcat(expwdir, 'EndAnaVals.dat' );

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
    
    if k == begi
       plot(x1,h1a,'--k','DisplayName', 'Analytic SWWE'); 
       hold on
    end
    
    plot(x1,h1,'-','DisplayName',strcat("\beta_2 = ",num2str(k*0.5) ));
end
legend('show');
hold off





