% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/REGSWWE/";
wdir = "AspRat1to0p01/DBalpha0p1/Beta";
dxdir = '/dx06/';

dx = 7.8126220722198776E-003;
scenbeg = -250;
scenend = 250;

linestarti = round(((-150 - scenbeg)/dx));
lineendi = round(((150 - scenbeg)/dx));
linesep = 1;

ksep = 1;
kend = 5;
figure;
hold on;
for k = 0:ksep:kend
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdirbase,wdir,BetaNumStr,dxdir);
    
    EndA = importdata(strcat(expdir,'EndAna.dat' ));
    End = importdata(strcat(expdir,'End.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));

    beta1str = extractBetween(param,"beta1 :","beta2");
    beta1 = str2double(beta1str{1,1});

    beta2str = extractAfter(param,"beta2 :");
    beta2 = str2double(beta2str);
    
    t = End(1,1);
    x = End(linestarti:linesep:lineendi,2);
    h = End(linestarti:linesep:lineendi,3);
    G = End(linestarti:linesep:lineendi,4);
    u = End(linestarti:linesep:lineendi,5);
    
    tA = EndA(linestarti:linesep:lineendi,1);
    xA = EndA(linestarti:linesep:lineendi,2);
    hA = EndA(linestarti:linesep:lineendi,3);
    GA = EndA(linestarti:linesep:lineendi,4);
    uA = EndA(linestarti:linesep:lineendi,5);
    
    plot(x,h,'-','DisplayName',strcat('\beta_1 =', num2str(beta1), '  \beta_2=',num2str(beta2)));
    xlabel('x (m)');
    ylabel('h (m)');
    
end
legend('show')
hold off


% EndA = importdata(strcat(wdir,'EndAna.dat' ));
% End = importdata(strcat(wdir,'End.dat' ));
% Init = importdata(strcat(wdir,'Init.dat' ));
% 
% t = End(:,1);
% x = End(:,2);
% h = End(:,3);
% G = End(:,4);
% u = End(:,5);
% 
% tA = EndA(:,1);
% xA = EndA(:,2);
% hA = EndA(:,3);
% GA = EndA(:,4);
% uA = EndA(:,5);
% 
% 
% figure;
% subplot(1,2,1);
% plot(x,h,'-b',x,hA,'--k');
% xlabel('x (m)');
% ylabel('h (m)');
% 
% subplot(1,2,2);
% plot(x,h - hA,'-r');
% xlabel('x (m)');
% ylabel('h - hA (m)');
% 
% figure;
% subplot(1,2,1);
% plot(x,u,'-b',x,uA,'--k');
% xlabel('x (m)');
% ylabel('u (m/s)');
% 
% subplot(1,2,2);
% plot(x,u - uA,'-r');
% xlabel('x (m)');
% ylabel('u - uA (m/s)');
% 
% 
% 
% figure;
% subplot(1,2,1);
% plot(x,G,'-b',x,GA,'--k');
% xlabel('x (m)');
% ylabel('G (m/s)');
% 
% subplot(1,2,2);
% plot(x,G - GA,'-r');
% xlabel('x (m)');
% ylabel('G - GA (m)');
