% Plot the results from gSerre

clc;
clear all;
close all;

load gSerre.r;
n = max(size(gSerre))/2;
h = gSerre(:,2);
uh = gSerre(:,3);
x = gSerre(:,1);

load gSerre.t
t = gSerre

figure(1)
plot(x(n+1:2*n),h(n+1:2*n),'or');
hold on
%axis([x(1) x(n) 0 ceil(max(h))]);
%axis([500 600 4 6]);
axis([x(1) x(n) 0.8 2]);

figure(2)
plot(x(n+1:2*n),uh(n+1:2*n),'or');
hold on
axis([x(1) x(n) 0 ceil(max(uh))]);

figure(3)
u = uh*0;
i = find(h~=0);
u(i) = uh(i)./h(i);
plot(x(n+1:2*n),u(n+1:2*n),'or');
hold on
axis([x(1) x(n) 0 ceil(max(u))]);


load gSerre.obs;
xobs = gSerre(:,1);
hobs = gSerre(:,2);
uobs = gSerre(:,3);


% plot flow depth
figure(1)
plot(xobs,hobs,'-b')
%axis([x(1) x(n) 0 ceil(max(h3))+1]);
xlabel('x')
ylabel('water depth')
title('Dam break problem, second-order')

%
% plot momentum
figure(2)
plot(xobs,uobs.*hobs,'-b')
%axis([x(1) x(n) 0 ceil(max(uh3))+5]);
xlabel('x')
ylabel('momentum')
title('Dam break problem, second-order')

%
% plot velocity
figure(3)
plot(xobs,uobs,'-b')
%axis([x(1) x(n) 0 ceil(max(u3))+1]);
xlabel('x')
ylabel('water velocity')
title('Dam break problem, second-order')

% %
% % Plot TVD results
% figure(4)
% load gSerre.TVD;
% n = max(size(gSerre));
% TVD = gSerre(:,2);
% x = gSerre(:,1);
% semilogx(x,TVD,'^r');







