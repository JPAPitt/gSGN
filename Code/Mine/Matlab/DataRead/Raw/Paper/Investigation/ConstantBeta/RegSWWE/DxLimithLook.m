% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/REGSWWE/";
wdir = "AspRat10to1/DBalpha0p5/Beta00";
dxdir = '/dx';


ksep = 1;
kend = 6;
figure;
hold on;
for k = 0:ksep:kend
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdirbase,wdir,dxdir,BetaNumStr,'/');
   
    End = importdata(strcat(expdir,'End.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));

    beta1str = extractBetween(param,"beta1 :","beta2");
    beta1 = str2double(beta1str{1,1});

    beta2str = extractAfter(param,"beta2 :");
    beta2 = str2double(beta2str);
    
    dxstr =  extractBetween(param," dx = (x_end - x_start) / (x_len -1)  :", "tstart");
    dx = str2double(dxstr{1,1});
    
    t = End(1,1);
    x = End(:,2);
    h = End(:,3);
    G = End(:,4);
    u = End(:,5);

    plot(x,h,'-','DisplayName',strcat('dx=',num2str(dx)));
%     plot(x,G,'-','DisplayName',strcat('dx=',num2str(dx)));

    
end

EndA = importdata(strcat(expdir,'EndAna.dat' ));

xA = EndA(:,2);
hA = EndA(:,3);
GA = EndA(:,4);
uA = EndA(:,5);

plot(xA,hA,'--k','DisplayName','SWWE Dambreak Solution');
% plot(xA,GA,'--k','DisplayName','SWWE Dambreak Solution');
xlabel('x (m)');
ylabel('h (m)');
% ylabel('G (m^2/s)');
title(strcat('\beta_2=',num2str(beta2),'  \beta_1 = -2/3 + \beta_2'))
legend('show')
hold off

%matlab2tikz('RegSWWEh.tex');
