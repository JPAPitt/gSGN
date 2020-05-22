% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/Reg2/SmoothDB/alpha2/timeseries/';


figure;

linesep = 20;
begi = 00;
sep = 2;
bege = 10;
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
    
    plot(x1(1:linesep:end),h1(1:linesep:end),'-','DisplayName',strcat("\beta_1 = ",num2str(beta1val,2), "  \beta_2 = ",num2str(beta2val,2)));
    
    if k == begi
       hold on
    end
end
hold off

legend('hide');
xlabel('x (m)');
ylabel('h (m)');
axis([-150 150 0.8 2.4]);
xticks([-150,-100,-50,0,50,100,150]);
yticks([0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4]);

matlab2tikz('SDBa2Reg2h.tex');

clc;
close all;

figure;
linesep = 20;
begi = 00;
sep = 2;
bege = 10;
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
    
    plot(x1(1:linesep:end),h1(1:linesep:end),'-','DisplayName',strcat("\beta_1 = ",num2str(beta1val,2), "  \beta_2 = ",num2str(beta2val,2)));
    
    if k == begi
       hold on
    end
end
hold off

legend('hide');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-150 150 -1 3]);
xticks([-150,-100,-50,0,50,100,150]);
yticks([-1,0,1,2,3]);

matlab2tikz('SDBa2Reg2u.tex');







