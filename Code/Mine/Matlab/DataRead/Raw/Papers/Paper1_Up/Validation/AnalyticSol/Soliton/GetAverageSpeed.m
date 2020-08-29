% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/SerreSoliton/";

a0 = 1.0;
a1 = 0.7;
g = 9.81;
c = sqrt(g*(a0 + a1));

n = 12;
SpeedUB = zeros(n+1,1);
SpeedLB = zeros(n+1,1);
dxs = zeros(n+1,1);

for k = 0:n
    dxNumStr = compose("%2.2d",k);
    expdir = strcat(wdir,dxNumStr,'/');
    EndF = importdata(strcat(expdir, 'End.dat'));
    t = EndF(1,1);
    x = EndF(:,2);
    h = EndF(:,3);
    G = EndF(:,4);
    u = EndF(:,5);
    
    param = fileread(strcat(expdir ,'Params.dat' ));
    dxstr = extractBetween(param,"dx","tstart");
    dx = str2double(dxstr{1,1});

    PeakLocj = find(h==max(h));
    
    xUB = x(PeakLocj) + 0.5*dx;
    xLB = x(PeakLocj) - 0.5*dx;
    
    SpeedUB(k+1) = xUB / t;
    SpeedLB(k+1) = xLB / t;
    dxs(k+1) = dx;


end

semilogx(dxs,SpeedUB,'x b',dxs,SpeedLB,'x r');
hold on;
plot([10^(-4),1],[c,c],'--k');
legend('hide');
xlabel('\Delta x');
ylabel('Peak Speed');
axis([10^-3 10 3.1 4.1]);
xticks([10^-3,10^-2,10^-1,10^0,10]);
yticks([3.1,3.3,3.5,3.7,3.9,4.1]);
matlab2tikz('PeakSpeedEstimates.tex');

