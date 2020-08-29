% Process Fortran Outputs

clc;
clear all;
close all;
format long g;

% Get list of directories to loop over when reading data
wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/SerreImpDispRedo/";
wdir = "AspRat2to1/DBalpha0p1/Beta";
dxdir = '/dx06/';

g= 9.81;
hl = 2;
hr = 1;

kstart = 0
ksep = 1;
kend = 2;
RelErrors = [];
for k = kstart:ksep:kend
    BetaNumStr = compose("%2.2d",k);
    expdir = strcat(wdirbase,wdir,BetaNumStr,dxdir);
   
    Energy = importdata(strcat(expdir,'Energy.dat' ));
    param = fileread(strcat(expdir ,'Params.dat' ));

    beta1str = extractBetween(param,"beta1 :","beta2");
    beta1 = str2double(beta1str{1,1});

    beta2str = extractAfter(param,"beta2 :");
    beta2 = str2double(beta2str);
    
    EndTimeStr = extractBetween(param,"actual end time :","dt");
    EndTime = str2double(EndTimeStr{1,1});
    
    InitialStr = extractAfter(Energy,"Initial");
    InitialEnergies = sscanf(InitialStr{2,1},'%f');
    
    EndStr = extractAfter(Energy,"End");
    SplitEndStr = split(EndStr{3,1});
    SplitEndStr = SplitEndStr(2:5);
    EndEnergies = str2double(SplitEndStr);
    
    uhflux = 0.5*g*(hl*hl - hr*hr)*EndTime;
    
    FluxEnergies = [0;uhflux;uhflux;0];
    
    InitialPlusFlux = InitialEnergies + FluxEnergies ;
    
    Error = EndEnergies - (InitialPlusFlux);
    RelError = Error ./ (InitialPlusFlux);
    RelError = [beta1;beta2;RelError];
    RelErrors = [RelErrors ,RelError ];

end

RelErrors = RelErrors';
