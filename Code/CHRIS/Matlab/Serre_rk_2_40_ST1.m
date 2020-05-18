% Plot the results from Serre_rk_2_40_Surface_Tension

clc;
clear all;
close all;

load Serre_rk_2_40_Surface_Tension.r;
n = max(size(Serre_rk_2_40_Surface_Tension))/2;
h = Serre_rk_2_40_Surface_Tension(:,2);
uh = Serre_rk_2_40_Surface_Tension(:,3);
x = Serre_rk_2_40_Surface_Tension(:,1);

load Serre_rk_2_40_Surface_Tension.t
t = Serre_rk_2_40_Surface_Tension

figure(1)
plot(x(n+1:2*n),h(n+1:2*n),'or');
hold on
%axis([x(1) x(n) 0 ceil(max(h))]);
%axis([500 600 4 6]);
%axis([x(1) x(n) 1.0 2]);

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


load Serre_rk_2_40_Surface_Tension.obs;
xobs = Serre_rk_2_40_Surface_Tension(:,1);
hobs = Serre_rk_2_40_Surface_Tension(:,2);
uobs = Serre_rk_2_40_Surface_Tension(:,3);


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
% load Serre_rk_2_40_Surface_Tension.TVD;
% n = max(size(Serre_rk_2_40_Surface_Tension));
% TVD = Serre_rk_2_40_Surface_Tension(:,2);
% x = Serre_rk_2_40_Surface_Tension(:,1);
% semilogx(x,TVD,'^r');







