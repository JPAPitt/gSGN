% Process Fortran Outputs

clc;
clear all;
close all;
format long g;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/Validation/AnalyticSolutions/DBSWWE/"

g= 9.81;
hl = 2;
hr = 1;

kstart = 0
ksep = 1;
kend = 7;
RelErrors = [];
for k = kstart:ksep:kend
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdir,BetaNumStr,'/');
   
    Energy = importdata(strcat(expdir,'Energy.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));

    Str1 = extractBetween(param,"dx","tstart ");
    dx = str2double(Str1{1,1});
    
    Str1 = extractBetween(param,"actual end time","dt");
    EndTime = str2double(Str1{1,1});
    
    InitialStr = extractAfter(Energy,"Initial");
    InitialEnergies = sscanf(InitialStr{2,1},'%f');
    
    EndStr = extractAfter(Energy,"End");
    SplitEndStr = split(EndStr{3,1});
    SplitEndStr = SplitEndStr(2:5);
    EndEnergies = str2double(SplitEndStr);
    
    uhflux = 0.5*g*(hl*hl - hr*hr)*EndTime;
    
    FluxEnergies = [0;uhflux;uhflux;3];
    
    InitialPlusFlux = InitialEnergies + FluxEnergies ;
    
    Error = EndEnergies - (InitialPlusFlux);
    RelError = abs(Error ./ (InitialPlusFlux));
    RelError = [dx;RelError];
    RelErrors = [RelErrors ,RelError ];

end

RelErrors = RelErrors';

dxs = RelErrors(:,1);
hs = RelErrors(:,2);
Gs = RelErrors(:,3);
uhs = RelErrors(:,4);
Hs = RelErrors(:,5);

figure;
loglog(dxs,hs,'s b',dxs,Gs,'o r',dxs,uhs,'^ k',dxs,Hs,'x g')
grid off
legend('h','G', 'uh', 'H','Location','northwest')
xlabel('\Delta x')
ylabel('C_1')
axis([10^-3 1 10^-16 1]);
xticks([10^-4,10^-3,10^-2,10^-1,10^0,10]);
yticks([10^-16,10^-12,10^-8,10^-4,1]);
matlab2tikz('EnergyResults.tex');
