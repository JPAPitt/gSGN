% Process Fortran Outputs

clc;
clear all;
close all;


kstart = 0;
ksep = 1;
kend = 3;
figure;
hold on;
for k = kstart:ksep:kend
    NumStr = compose("%2.2d",k);
    TimeSeries = importdata(strcat('Soliton_Ex',NumStr,'.dat'));
    t = TimeSeries(:,1);
    H = TimeSeries(:,6);
    dx= TimeSeries(1,2);
    plot(t,H,'.');
    
end

xlabel('t (s)');
ylabel('H (m^3/s^2)');
legend('Serre', 'Serre Improved Dispersion', 'SWWE', 'RegSWWE');
matlab2tikz('EnergyOverTime.tex');


    
