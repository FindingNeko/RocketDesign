%MAE 5391
%Bryan Rathke
clc;clear all;close all
%Constants
c=7000;
mdot=280;
Qr=2400;
jay=777.9;
Pchem=mdot*Qr*jay;
v=linspace(0,21000,21001);
vrat=v./c;
Pjet=(1/2)*mdot.*c.^2;

%Propulsive efficiency
etaP=(2.*(vrat))./(1.+(vrat.^2));

%Internal efficiency
etaI=Pjet/Pchem;
eatI=linspace(1,1,21001).*etaI;

%Overall efficiency
etaOA=etaP.*etaI;


figure(1)
subplot(3,1,1)
plot(vrat,etaP)
title('Propulsive efficiency')
subplot(3,1,2)
plot(vrat,etaI)
title('Internal efficiency')
subplot(3,1,3)
plot(vrat,etaOA)
title('Overall efficiency')