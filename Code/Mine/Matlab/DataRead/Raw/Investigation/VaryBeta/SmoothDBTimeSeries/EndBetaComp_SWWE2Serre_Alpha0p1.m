% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/SWWE2Serre/SmoothDB/alpha0p1/timeseries/';


figure;

linesep = 20;

begi = 00;
sep =1;
bege = 3;
for k = begi:sep:bege
    if k < 10
        expwdir = strcat(wdir,'0',int2str(k), '/');
    else 
        expwdir = strcat(wdir,int2str(k), '/');
    end 

    endfile = strcat(expwdir, 'EndVals.dat' );
    endanafile = strcat(expwdir, 'EndAnaVals.dat' );
    paramfile = strcat(expwdir, 'Params.dat' );

    xhGuend = importdata(endfile);
    xhGuendana = importdata(endanafile);
    param = fileread(paramfile);
    
    beta1 = extractBetween(param,"beta1 : ","beta2");
    beta1val = str2double(beta1{1,1});
    beta2 = extractAfter(param,"beta2 :");
    beta2val = str2double(beta2);

    x1 = xhGuend(:,1);
    h1 = xhGuend(:,2);

    x1a = xhGuendana(:,1);
    h1a = xhGuendana(:,2);
    
    plot(x1(1:linesep:end),h1(1:linesep:end),'-','DisplayName',strcat("\beta_1 =",num2str(beta1val,2)));
    
    if k == begi
       hold on
    end
end
legend('show');
xlabel('x (m)');
ylabel('h (m)');
axis([-100 100 0.9 2.1]);
xticks([-100,-50,0,50,100]);
yticks([1,1.25,1.5,1.75,2]);
hold off
matlab2tikz('SDBBetaAroundSWWEh.tex');

clc;
close all;


linesep = 20;

figure;
begi = 00;
sep = 1;
bege = 3;
for k = begi:sep:bege
    if k < 10
        expwdir = strcat(wdir,'0',int2str(k), '/');
    else 
        expwdir = strcat(wdir,int2str(k), '/');
    end 

    endfile = strcat(expwdir, 'EndVals.dat' );
    endanafile = strcat(expwdir, 'EndAnaVals.dat' );
    paramfile = strcat(expwdir, 'Params.dat' );

    xhGuend = importdata(endfile);
    xhGuendana = importdata(endanafile);
    param = fileread(paramfile);
    
    beta1 = extractBetween(param,"beta1 : ","beta2");
    beta1val = str2double(beta1{1,1});
    beta2 = extractAfter(param,"beta2 :");
    beta2val = str2double(beta2);

    x1 = xhGuend(:,1);
    h1 = xhGuend(:,4);

    x1a = xhGuendana(:,1);
    h1a = xhGuendana(:,4);
    
    plot(x1(1:linesep:end),h1(1:linesep:end),'-','DisplayName',strcat("\beta_1 =",num2str(beta1val,2)));
    
    if k == begi
       hold on
    end
end
legend('show');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-100 100 0 2.5]);
xticks([-100,-50,0,50,100]);
yticks([0,0.5,1,1.5,2,2.5]);
hold off
matlab2tikz('SDBBetaAroundSWWEu.tex');






