% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolDBLoop/theta1/08/";
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/SWWE/AnaSolDBLoop/theta1/07/"

wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/";
wdir = "/DXLook/ImpDispSerre/AspRat2to1/DBalpha0p1/dx";

dx = 7.8126220722198776E-003;
scenbeg = -180;
scenend = 180;


kstart = 5;
ksep = 1;
kend = 8;
figure;
hold on;
for k = kstart:ksep:kend
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdirbase,wdir,BetaNumStr,'/');
   
    End = importdata(strcat(expdir,'End.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));

    dxstr = extractBetween(param,"dx","tstart");
    dx = str2double(dxstr{1,1});
    
    t = End(1,1);
    x = End(:,2);
    h = End(:,3);
    G = End(:,4);
    u = End(:,5);

    plot(x,u,'-','DisplayName',strcat('\beta_2=',num2str(dx)));
    xlabel('x (m)');
    ylabel('h (m)');
    
end

EndA = importdata(strcat(expdir,'EndAna.dat' ));
tA = EndA(1,1);
xA = EndA(:,2);
hA = EndA(:,3);
GA = EndA(:,4);
uA = EndA(:,5);

plot(xA,uA,'--k','DisplayName','SWWE Dambreak Solution');
legend('show');

% axis([-175 175 1 2])
% xticks([-175,-87.5,0,87.5,175])
% yticks([1,1.2,1.4,1.6,1.8,2])
% legend('show')
% hold off


% matlab2tikz('h.tex');
