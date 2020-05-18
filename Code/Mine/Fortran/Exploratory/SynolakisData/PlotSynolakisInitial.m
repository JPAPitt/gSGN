% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wfile = 'CellNodesAll.dat';

xhwub = importdata(wfile);

x = xhwub(:,1);
h = xhwub(:,2);
w = xhwub(:,3);
u = xhwub(:,4);
b = xhwub(:,5);

f1 = figure;
sgtitle('All Quantities');
subplot(2, 2, 1);
plot(x,h,'-g');
title('h') ;

subplot(2,2, 2);
plot(x,w,'-b');
title('w') ;

subplot(2, 2, 3);
plot(x,u,'-r');
title('u') ;

subplot(2, 2, 4);
plot(x,b,'-k');
title('b') ;

f2 = figure;
sgtitle('Water and Velocity Profiles');
subplot(1, 2, 1);
plot(x,w,'-b',x,b,'-k');
title('water profile') ;
legend('w','b');
subplot(1, 2, 2);
plot(x,u,'-r');
title('velocity profile');

f3 = figure;
sgtitle('Zoom In Profiles');
subplot(1, 2, 1);
plot(x,w,'-b',x,b,'-k');
xlim([-10 70]);
ylim([-0.05,0.05]);
title('water profile') 
legend('w','b');
subplot(1, 2, 2);
plot(x,u,'-r');
xlim([-10 70]);
title('velocity profile'); 

