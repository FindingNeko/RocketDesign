%Liquid_Rocket_Design
%Bryan Rathke
clear all, close all
%% Historical data

% Dim_v=linspace(0,60);
% F_v=linspace(0,4000);
% Length_v=0.005048*F_v+31.92;
% D_v=0.00357*F_v+14.48;
% 
% figure(1)
% hold on
% plot(F_v,Length_v)
% plot(F_v,D_v)
% xlabel('Thrust (F)')
% ylabel('Dimension (cm)')

%% Design parameters

DeltaV=1721;        %Required change in velocity [m/s]
Mpay=4914;          %Payload mass [kg]
TTWR=0.3;           %Thrust to weight ratio (minimum)
Mnought=12000;      %Initial vehicle mass [kg]
EnvH=300;           %Design envelope hieght [cm]
EnvW=300;           %Design envelope width [cm]
Pc=700000;          %Chamber pressure [Pa]
Pa=0;               %Atmosphere
Ncstar=1;           %CStar combustion efficiency
lambda=0.98;        %Nozzle Efficiency

g=9.80665;          %Standard gravity

%Combustion parameters

OTFR=2.3;
FlameT=3510;   
k=1.225;
MM=22.1;
R=8314/MM;          %Gas constant
rho_f=810;
rho_o=1142;

%% Initial calculations, space bipropellant

DeltaV0=DeltaV*1.1;
F=Mnought*g*TTWR;
Meng=F/(g*(0.000609*F+13.44));
Leng=0.0054*F+31.92;
Deng=0.00357*F+14.48;

if (Leng>EnvH) || (Deng>EnvW);          %Check engine size vs envelope
    print('Warning: Estimated engine size exceeds design envelope.')
end

kmo=k-1;
kpo=k+1;
omk=1-k;

cstar=(Ncstar*(k*R*FlameT)^0.5) / ...
    (k*(2/kpo)^(kpo/(2*kmo)));

Pe_v=linspace(0,45000,1000);                           %Exit pressure vector
Me_v=(((Pe_v/Pc).^(omk/k)-1)/((k-1)/2)).^(0.5);  %Exit Mach Number
epsilon_v=(Me_v.^-1).*((2/kpo).*(1+(kmo/2).*Me_v.^2)).^(kpo/(2*kmo));
Isp_v=lambda.*(((cstar*k)/g).*((2/kmo).*(2/kpo)...
    .^(kpo/kmo).*(1-(Pe_v./Pc).^(kmo/k))).^0.5 ... 
    +((cstar.*epsilon_v)./(g*Pc)).*(Pe_v-Pa));
figure(1)
plot(epsilon_v,Isp_v)
grid on
title('Space engine performance')
xlabel('Nozzle Expansion Ratio')
ylabel('Specific Impulse Isp')

prompt1='Select a nozzle expansion ratio: ';
epsilon=input(prompt1);
if isempty(epsilon)
    epsilon=100
end

%Isp=lambda*(((Cstar*k)/g)*((2/kmo)*(2/kpo)...
%    ^(kpo/kmo)*(1-(Pe_v/Pc)^(kmo/k)))^0.5 ... 
%    +((Cstar.*epsilon_v)./(g*Pc)).*(Pe_v-Pa));

%% Pressure calculations

V_fuel=10;
V_ox=10;
PD_inj=Pc*0.2;
PD_feed=50000;
PD_dyn_f=0.5*rho_f*V_fuel^2;
PD_dyn_o=0.5*rho_o*V_ox^2;

P_ftank=Pc+PD_inj+PD_feed+PD_dyn_f;
P_otank=Pc+PD_inj+PD_feed+PD_dyn_o;



%Pt_pres=(10.^(-0.1068.*(log(Vt)-0.2588))).*(10^6);

%% Fuel mass


Isp=337;                    %Assumed
Mtotal=Mnought+1;
c=.5;
while(Mtotal>Mnought)
    f_inert=c;
    Mpro=(Mpay*(exp(DeltaV0/(Isp*g))-1)*(1-f_inert))/(1-(f_inert*exp(DeltaV0/(Isp*g))));
    Mtotal=(Mpro/(1-f_inert))+Mpay;
    c=c-0.001;
end
Mpro
Mtotal
Mo=(OTFR*Mpro/(OTFR+1));
Mf=Mpro/(1+OTFR);

%% Tank mass

Vt_ox=1.1*(Mo/rho_o);
Vt_f=1.1*(Mf/rho_f);
mdot=F/(Isp*g);
At=mdot*cstar/Pc;
Ae=epsilon*At;
Dt=(4*At/pi)^.5;
De=(4*Ae/pi)^.5;

prompt2='Select a characteristic length for C.C. : ';
Lstar=input(prompt2);
if isempty(Lstar)
    Lstar=0.9
end
prompt3='Select a chamber mach number: ';
Machc=input(prompt3);
if isempty(Machc)
    Machc=0.2
end
Ac=(At/Machc)*((2/kpo)*(1+(kmo/2)*Machc^2))^(kpo/(2*kmo));
Lc=Lstar*At/Ac;
Contraction_Ratio=Ac/At
Dc=(4*Ac/pi)^.5
prompt4='Select a multiplication factor: ';
Phic=input(prompt4);
if isempty(Phic)
    Phic=3
end
prompt5='Enter yield strength of C.C. material: ';
prompt6='Enter CC material density: ';;
Ftu=input(prompt5);
if isempty(Ftu)
    Material='Columbium'
    Ftu=310000000
    rho_c=8500
else
    rho_c=input(prompt)
end
prompt7='Enter contraction angle in degrees: ';
Thetac=input(prompt7);
if isempty(Thetac)
    Thetac=45
end
cwall=(Pc*Dc*Phic)/(2*Ftu);
Mc=pi*rho_c*cwall*(Dc*Lc+(pi*((Dc^2-Dt^2)/4)/tand(Thetac)));
prompt8='Choose a nozzle cone half-angle: ';
Thetan=input(prompt8);
if isempty(Thetan)
    Thetan=15
end
Ln=(De-Dt)/(2*tand(Thetan));
prompt9='Enter "percent bell" fraction of 15 degree conical nozzle: '
if (epsilon>29)
    if lambda==0.99
        Lf=0.86;
    elseif lambda==.987
        Lf=0.8;
    elseif lambda==.985
        Lf=0.75;
    elseif lambda==.98
        Lf=0.675;
    else
        Lf=input(prompt9);
    end
else
    Lf=input(prompt9);
end

%% Nozzle drawing
Rt=Dt/2;
Re=De/2;
prompt10='Enter initial parabola angle: ';
prompt11='Enter final parabola angle: ';
Thetaip=input(prompt10);
Thetafp=input(prompt11);
% x0=linspace((-1.5*Rt),0);
% x1=linspace(0,(Rt));
% x2=linspace((1.376*Rt),Ln);
% Noz0=(2.5*Rt)-(((1.5*Rt)^2)-x0.^2);
% Noz1=(2.382*Rt)-(((.882*Rt)^2)-x1.^2);
% %Noz2
% figure(2)
% plot(x0,Noz0)
% hold on;
% plot(x1,Noz1)
% xlim([-2*Rt,(1.5*Ln)])
% ylim([0,(1.5*Re)])
figure(2)
function h = circle(x,y,r,th)
hold on
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end
circle
