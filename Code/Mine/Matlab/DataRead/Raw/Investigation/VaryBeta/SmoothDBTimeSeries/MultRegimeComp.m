% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/Models/gSGN/ConstantBetas/Serre/Soliton/';


wdirSerre = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/ImpDisp/SmoothDB/alpha2/timeseries/00/';

initfile = strcat(wdirSerre, 'InitVal.dat' );

%make initial condition figure
figure;
xhGuinit = importdata(initfile);
x0 = xhGuinit(:,1);
h0 = xhGuinit(:,2);
G0 = xhGuinit(:,3);
u0 = xhGuinit(:,4);
subplot(2, 2, 1);
plot(x0,h0);
xlabel('x (m)');
ylabel('h (m)');
axis([-100 100 0.9 2.1]);

subplot(2, 2, 2);
plot(x0,h0);
xlabel('x (m)');
ylabel('h (m)');
axis([-5 5 0.9 2.1]);


subplot(2, 2, 3);
plot(x0,u0);
xlabel('x (m)');
ylabel('u (m/s)');
axis([-100 100 -0.1 0.1]);

subplot(2, 2, 4);
plot(x0,G0);
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-100 100 -0.1 0.1]);

sgtitle('Initial Conditions - Smooth Dambreak \alpha = 0.1');


wdirSWWE = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/RegSWWE/SmoothDB/alpha2/timeseries/00/';
wdirRegSWWE = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimhG/VaryBeta/RegSWWE/SmoothDB/alpha2/timeseries/03/';

wdirs = {wdirSerre,wdirSWWE,wdirRegSWWE};

figure;
subplot(1,4,1)
for k = 1:3
    endfile = strcat(wdirs{1,k}, 'EndVals.dat' );
    paramfile = strcat(wdirs{1,k}, 'Params.dat' );
    
    xhGuend = importdata(endfile);
    param = fileread(paramfile);
    
    beta1 = extractBetween(param,"beta1 : ","beta2");
    beta1val = str2double(beta1{1,1});
    beta2 = extractAfter(param,"beta2 :");
    beta2val = str2double(beta2);

    x1 = xhGuend(:,1);
    h1 = xhGuend(:,2);
    
    plot(x1,h1,'-', 'DisplayName',strcat("\beta_1 = ",rat(beta1val) , ' and \beta_2 = ',rat(beta2val) ));
    
    if k == 1
       hold on
    end
    
end 
xlabel('x (m)');
ylabel('h (m)');
axis([-100 100 0.9 2.1]);
legend('show');
hold off


subplot(1,4,2)
for k = 1:3
    endfile = strcat(wdirs{1,k}, 'EndVals.dat' );
    paramfile = strcat(wdirs{1,k}, 'Params.dat' );
    
    xhGuend = importdata(endfile);
    param = fileread(paramfile);
    
    beta1 = extractBetween(param,"beta1 : ","beta2");
    beta1val = str2double(beta1{1,1});
    beta2 = extractAfter(param,"beta2 :");
    beta2val = str2double(beta2);

    x1 = xhGuend(:,1);
    G1 = xhGuend(:,3);
    
    plot(x1,G1,'-', 'DisplayName',strcat("\beta_1 = ",rat(beta1val) , ' and \beta_2 = ',rat(beta2val) ));
    
    if k == 1
       hold on
    end
    
end 
xlabel('x (m)');
ylabel('G (m^2/s)');
%axis([-100 100 0.9 2.1]);
legend('show');
hold off

subplot(1,4,3)
for k = 1:3
    endfile = strcat(wdirs{1,k}, 'EndVals.dat' );
    paramfile = strcat(wdirs{1,k}, 'Params.dat' );
    
    xhGuend = importdata(endfile);
    param = fileread(paramfile);
    
    beta1 = extractBetween(param,"beta1 : ","beta2");
    beta1val = str2double(beta1{1,1});
    beta2 = extractAfter(param,"beta2 :");
    beta2val = str2double(beta2);

    x1 = xhGuend(:,1);
    u1 = xhGuend(:,4);
    
    plot(x1,u1,'-', 'DisplayName',strcat("\beta_1 = ",rat(beta1val) , ' and \beta_2 = ',rat(beta2val) ));
    
    if k == 1
       hold on
    end
    
end 
xlabel('x (m)');
ylabel('u (m/s)');
%axis([-100 100 0.9 2.1]);
legend('show');
hold off


subplot(1,4,4)
for k = 1:3
    endfile = strcat(wdirs{1,k}, 'EndVals.dat' );
    paramfile = strcat(wdirs{1,k}, 'Params.dat' );
    
    xhGuend = importdata(endfile);
    param = fileread(paramfile);
    
    beta1 = extractBetween(param,"beta1 : ","beta2");
    beta1val = str2double(beta1{1,1});
    beta2 = extractAfter(param,"beta2 :");
    beta2val = str2double(beta2);

    x1 = xhGuend(:,1);
    h1 = xhGuend(:,2);
    u1 = xhGuend(:,4);
    
    plot(x1,u1 .* h1,'-', 'DisplayName',strcat("\beta_1 = ",rat(beta1val) , ' and \beta_2 = ',rat(beta2val) ));
    
    if k == 1
       hold on
    end
    
end 
xlabel('x (m)');
ylabel('uh (m^2/s)');
%axis([-100 100 0.9 2.1]);
legend('show');
hold off





