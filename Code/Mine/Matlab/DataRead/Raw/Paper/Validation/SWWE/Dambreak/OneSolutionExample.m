% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/Validation/AnalyticSolutions/DBSWWE/04/"

param = fileread(strcat(wdir ,'Params.dat' ));
dxstr = extractBetween(param,"dx","tstart");
dx = str2double(dxstr{1,1});

scenbeg = -250;
scenend = 250;


lessdenseareab = -175;
denseareab =145;
denseareae =150;
lessdenseareae = 175;

denselinesep = 1;
densedotsep = 1;
lessdenselinesep = 50;
lessdensedotsep = 50;



lessdenseareabi = round(((lessdenseareab - scenbeg)/dx));
denseareabi = round(((denseareab - scenbeg)/dx));
denseareaei = round(((denseareae - scenbeg)/dx));
lessdenseareaei = round(((lessdenseareae - scenbeg)/dx));


EndA = importdata(strcat(wdir,'EndAna.dat' ));
End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));


t0 = Init(1,1);
x0 = [Init(lessdenseareabi:lessdenselinesep:denseareabi,2); ...
    Init(denseareabi:denselinesep:denseareaei,2); ...
    Init(denseareaei:lessdenselinesep:lessdenseareaei,2)];
h0 = [Init(lessdenseareabi:lessdenselinesep:denseareabi,3); ...
    Init(denseareabi:denselinesep:denseareaei,3); ...
    Init(denseareaei:lessdenselinesep:lessdenseareaei,3)];
G0 = [Init(lessdenseareabi:lessdenselinesep:denseareabi,4); ...
    Init(denseareabi:denselinesep:denseareaei,4); ...
    Init(denseareaei:lessdenselinesep:lessdenseareaei,4)];
u0 = [Init(lessdenseareabi:lessdenselinesep:denseareabi,5); ...
    Init(denseareabi:denselinesep:denseareaei,5); ...
    Init(denseareaei:lessdenselinesep:lessdenseareaei,5)];


t = End(1,1);
x = [End(lessdenseareabi:lessdensedotsep:denseareabi,2); ...
    End(denseareabi:densedotsep:denseareaei,2); ...
    End(denseareaei:lessdensedotsep:lessdenseareaei,2)];
h = [End(lessdenseareabi:lessdensedotsep:denseareabi,3); ...
    End(denseareabi:densedotsep:denseareaei,3); ...
    End(denseareaei:lessdensedotsep:lessdenseareaei,3)];
G = [End(lessdenseareabi:lessdensedotsep:denseareabi,4); ...
    End(denseareabi:densedotsep:denseareaei,4); ...
    End(denseareaei:lessdensedotsep:lessdenseareaei,4)];
u = [End(lessdenseareabi:lessdensedotsep:denseareabi,5); ...
    End(denseareabi:densedotsep:denseareaei,5); ...
    End(denseareaei:lessdensedotsep:lessdenseareaei,5)];

tA = EndA(1,1);
xA = [EndA(lessdenseareabi:lessdenselinesep:denseareabi,2); ...
    EndA(denseareabi:denselinesep:denseareaei,2); ...
    EndA(denseareaei:lessdenselinesep:lessdenseareaei,2)];
hA = [EndA(lessdenseareabi:lessdenselinesep:denseareabi,3); ...
    EndA(denseareabi:denselinesep:denseareaei,3); ...
    EndA(denseareaei:lessdenselinesep:lessdenseareaei,3)];
GA = [EndA(lessdenseareabi:lessdenselinesep:denseareabi,4); ...
    EndA(denseareabi:denselinesep:denseareaei,4); ...
    EndA(denseareaei:lessdenselinesep:lessdenseareaei,4)];
uA = [EndA(lessdenseareabi:lessdenselinesep:denseareabi,5); ...
    EndA(denseareabi:denselinesep:denseareaei,5); ...
    EndA(denseareaei:lessdenselinesep:lessdenseareaei,5)];



figure;
plot(x,h,'.r');
hold on;
plot(xA,hA,'-b');
plot(x0,h0,'-k');
hold off
xlabel('x (m)');
ylabel('h (m)');
axis([-175 175  0.9 2.1]);
xticks([-175,-87.5,0,87.5,175]);
yticks([0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1]);
legend('hide')
matlab2tikz('h.tex');
close all;


figure;
plot(x,G,'.r');
hold on;
plot(xA,GA,'-b');
plot(x0,G0,'-k');
hold off
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-175 175  0 2]);
xticks([-175,-87.5,0,87.5,175]);
yticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2]);
legend('hide')
matlab2tikz('G.tex');
close all;


figure;
plot(x,u,'.r');
hold on;
plot(xA,uA,'-b');
plot(x0,u0,'-k');
hold off
xlabel('x (m)');
ylabel('u (m/s)');
axis([-175 175  0 1.4]);
xticks([-175,-87.5,0,87.5,175]);
yticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4]);
legend('hide');
matlab2tikz('u.tex');
close all;

