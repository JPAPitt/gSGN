% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wfile = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/ReconFDM/Norms.dat';

Out = importdata(wfile);

dxs = Out(:,1);
Normhs = Out(:,2);
NormGs = Out(:,3);
Normus = Out(:,4);
Normdus = Out(:,5);


slope1 = 0.1*dxs.^(1);
slope2 = 0.1*dxs.^(2);
loglog(dxs,Normhs,'or',dxs,NormGs,'.g',dxs,Normus,'ob',dxs,Normdus,'.k' ,dxs,slope1,'-k',dxs,slope2,'-b')
grid on
legend('h','G', 'u','du','Slope 1', 'Slope 2','Location','northwest')



