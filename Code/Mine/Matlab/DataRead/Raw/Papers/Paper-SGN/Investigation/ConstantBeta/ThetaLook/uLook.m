% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/LimiterLook/";
wdir = "SWWE/AspRat2to1/DBalpha0p1/Theta";
dxdir = '/dx06/';

dx = 7.8126220722198776E-003;
scenbeg = -250;
scenend = 250;

linestart = -200;
lineend = 200;
linestarti = round(((linestart - scenbeg)/dx));
lineendi = round(((lineend - scenbeg)/dx));
linesep = 1;

kstart = 0;
ksep = 1;
kend = 20;
figure;
hold on;
for k = [kstart:ksep:kend]
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdirbase,wdir,BetaNumStr,dxdir);
   
    End = importdata(strcat(expdir,'End.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));
    
    thetastr = extractBetween(param,"theta","gravity");
    theta = str2double(thetastr{1,1});

    beta1str = extractBetween(param,"beta1","beta2");
    beta1 = str2double(beta1str{1,1});

    beta2str = extractAfter(param,"beta2");
    beta2 = str2double(beta2str);
    
    t = End(1,1);
    x = End(linestarti:linesep:lineendi,2);
    h = End(linestarti:linesep:lineendi,3);
    G = End(linestarti:linesep:lineendi,4);
    u = End(linestarti:linesep:lineendi,5);

    plot(x,u,'-','DisplayName',strcat('\theta=',num2str(theta)));
    
end

EndA = importdata(strcat(expdir,'EndAna.dat' ));
tA = EndA(1,1);
xA = EndA(linestarti:linesep:lineendi,2);
hA = EndA(linestarti:linesep:lineendi,3);
GA = EndA(linestarti:linesep:lineendi,4);
uA = EndA(linestarti:linesep:lineendi,5);

plot(xA,uA,'--k','DisplayName','SWWE Dambreak Solution');

xlabel('x (m)');
ylabel('u (m/s)');
legend('show')
hold off

% axis([135 155 0.0 7])
% xticks([135,140,145,150,155])
% yticks([0,1,2,3,4,5,6,7])
% legend('hide')
% hold off
% 
% matlab2tikz('GFront.tex');
