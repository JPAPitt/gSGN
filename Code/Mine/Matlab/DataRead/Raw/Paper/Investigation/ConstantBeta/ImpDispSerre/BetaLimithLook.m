% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/SerreImpDisp/";
wdir = "AspRat2to1/DBalpha0p1/Beta";
dxdir = '/dx06/';

dx = 7.8126220722198776E-003;
scenbeg = -250;
scenend = 250;

linestart = -175;
lineend = 175;
linestarti = round(((linestart - scenbeg)/dx));
lineendi = round(((lineend - scenbeg)/dx));
linesep = 1;

ksep = 5;
kend = 20;
figure;
hold on;
for k = 0:ksep:kend
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdirbase,wdir,BetaNumStr,dxdir);
   
    End = importdata(strcat(expdir,'End.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));

    beta1str = extractBetween(param,"beta1 :","beta2");
    beta1 = str2double(beta1str{1,1});

    beta2str = extractAfter(param,"beta2 :");
    beta2 = str2double(beta2str);
    
    t = End(1,1);
    x = End(:,2);
    h = End(:,3);
    G = End(:,4);
    u = End(:,5);

    plot(x,u,'-','DisplayName',strcat('\beta=',num2str(beta2)));
%     plot(x,G,'-','DisplayName',strcat('\beta_2=',num2str(beta2),'  \beta_1 = -2/3 + \beta_2'));

    
end

EndA = importdata(strcat(expdir,'EndAna.dat' ));

xA = EndA(:,2);
hA = EndA(:,3);
GA = EndA(:,4);
uA = EndA(:,5);

plot(xA,uA,'--k','DisplayName','SWWE Dambreak Solution');
% plot(xA,GA,'--k','DisplayName','SWWE Dambreak Solution');
xlabel('x (m)');
ylabel('h (m)');
% ylabel('G (m^2/s)');
legend('show')
hold off

%matlab2tikz('RegSWWEh.tex');
