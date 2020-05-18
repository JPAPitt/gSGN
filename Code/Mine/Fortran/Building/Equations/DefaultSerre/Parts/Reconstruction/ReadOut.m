% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wfile = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/ReconSine/Norms.dat';

Out = importdata(wfile);

dxs = Out(:,1);
Norms = Out(:,2);


slope1 = 0.1*dxs.^(1);
slope2 = 0.1*dxs.^(2);
loglog(dxs,Norms,'or',dxs,slope1,'-k',dxs,slope2,'-b')
grid on
legend('Norms','Slope 1', 'Slope 2','Location','northwest')



