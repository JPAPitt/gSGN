% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdirbase = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/BetaVary/SmoothDBInvestigation/";
wdirend = 'OneOffTests/10to1/dbalpha0p5/SWWE2Serre/';
wdir = strcat(wdirbase,wdirend);

times = [21,41,101,180];

for k = 1:size(times,2)
    timek = times(k);
    SimA = importdata(strcat(wdir, compose("%3d",timek) ,'Ana.dat' ));
    Sim = importdata(strcat(wdir, compose("%3d",timek) ,'.dat' ));
    param = fileread(strcat(wdir ,'Params.dat' ));

    t = Sim(1,1);
    x = Sim(:,2);
    h = Sim(:,3);
    G = Sim(:,4);
    u = Sim(:,5);
    xA = SimA(:,2);
    hA = SimA(:,3);
    GA = SimA(:,4);
    uA = SimA(:,5);

    b10str = extractBetween(param,"b10 :","b11");
    b10 = str2double(b10str{1,1});

    b11str = extractBetween(param,"b11 :","b20");
    b11 = str2double(b11str{1,1});

    b20str = extractBetween(param,"b20 :","b21");
    b20 = str2double(b20str{1,1});

    b21str = extractBetween(param,"b21 :","t0");
    b21 = str2double(b21str{1,1});

    t0str = extractAfter(param,"t0 :");
    t0 = str2double(t0str);

    figure;
    subplot(1,2,1);
    plot(x,h,'-b',x,hA,'--k');
    xlabel('x (m)');
    ylabel('h (m)');

    subplot(1,2,2);
    plot(x,G,'-r',x,GA,'--k');
    xlabel('x (m)');
    ylabel('G (m^2/s)');
    sgtitle(strcat('Initial \beta_1 = ',num2str(b10),' \beta_2 = ',num2str(b20) , '  switch time = ' ,num2str(t0),'  final \beta_1 = ',num2str(b11),' \beta_2 = ',num2str(b21),'  at t =',num2str(t)));
end
