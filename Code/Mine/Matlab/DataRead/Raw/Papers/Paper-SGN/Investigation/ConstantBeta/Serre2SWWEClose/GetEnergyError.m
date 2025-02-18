% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/SmoothDBInvestigation/Serre2SWWEClose/";
wdir = "AspRat2to1/DBalpha0p1/Beta";
dxdir = '/dx06/';

g= 9.81;
hl = 2;
hr = 1;

ksep = 10;
kend = 20;
RelErrors = [];
for k = 0:ksep:kend
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
    EndEnergies = sscanf(EndStr{3,1},'%f');
    
    uhflux = g/2*(hl^2 - hr^2)*EndTime;
    
    FluxEnergies = [0;uhflux;uhflux;0];
    
    Error = EndEnergies - (InitialEnergies +FluxEnergies );
    RelError = Error ./ ([InitialEnergies(1,1);1;1;InitialEnergies(4,1)]);
    RelError = [beta1;beta2;RelError];
    RelErrors = [RelErrors ,RelError ];

end

RelErrors = RelErrors';
