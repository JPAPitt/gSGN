% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/Energy/Soliton/";


kstart = 0;
ksep = 1;
kend = 3;
figure;
hold on;
for k = kstart:ksep:kend
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdir,BetaNumStr,'/');
   
    End = importdata(strcat(expdir,'End.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));

    beta1str = extractBetween(param,"beta1","beta2");
    beta1 = str2double(beta1str{1,1});

    beta2str = extractAfter(param,"beta2");
    beta2 = str2double(beta2str);
    
    t = End(1,1);
    x = End(:,2);
    h = End(:,3);
    G = End(:,4);
    u = End(:,5);

    plot(x,G,'-','DisplayName',strcat('\beta_1=',num2str(beta1),'\beta_2=',num2str(beta2)));
    xlabel('x (m)');
    ylabel('G (m^2/s)');
    
end
title(strcat('t = ',num2str(t)))
legend('show')
% 
% EndA = importdata(strcat(expdir,'EndAna.dat' ));
% tA = EndA(1,1);
% xA = EndA(:,2);
% hA = EndA(:,3);
% GA = EndA(:,4);
% uA = EndA(:,5);
% title(strcat('t = ',num2str(t)))
% plot(xA,hA,'--k','DisplayName','Soliton Serre');
% legend('show')

