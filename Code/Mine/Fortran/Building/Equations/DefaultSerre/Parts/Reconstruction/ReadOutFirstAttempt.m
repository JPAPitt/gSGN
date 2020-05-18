% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
DataDir = '/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/CodeAndData/Data/RAW/FortranTests/ReconSine';
DataDirs = dir(DataDir);
DataDirs = extractfield(DataDirs,'name');
shape = size(DataDirs);
DataDirs = DataDirs(3:shape(2));

Norms = zeros(size(DataDirs));
dxs = zeros(size(DataDirs));
for i = 1 : max(size(DataDirs))
    wdir = strcat(DataDir,'/',DataDirs{i},'/');
    wfile = strcat(wdir,'R.dat');
    Out = importdata(wfile);

    n = max(size(Out));
    x = Out(:,1);
    hA = Out(:,2);
    hN = Out(:,3);

    L2 = norm(abs(hA - hN),2) / norm(abs(hA),2);
    Norms(i) = L2;
    dxs(i) = x(6);
    
    clear n x hA hN L2 Out
end

slope1 = 0.1*dxs.^(1);
slope2 = 0.1*dxs.^(2);
loglog(dxs,Norms,'or',dxs,slope1,'-k',dxs,slope2,'-b')
grid on
legend('Norms','Slope 1', 'Slope 2','Location','northwest')

WriteOut = [dxs;Norms];
fileID = fopen('Norms.txt','w');
fprintf(fileID, 'Norms\n\n');
fprintf(fileID,'%e %e\n',WriteOut);


