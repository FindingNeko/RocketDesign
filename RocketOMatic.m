%Liquid_Space_Rocket_Design
%Bryan Rathke
clear all, close all

%% Design parameters
prompt0='Press enter at prompts to use default values: (press enter)';
defaults=input(prompt0);

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

OTFR=2.3;           %Combustion parameters
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

Tburn=Isp
%% CC and Nozzle 
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
prompt9='Enter "percent bell" fraction of 15 degree conical nozzle: ';
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
Lb=Ln*Lf;
%% Nozzle drawing
Rt=Dt/2;
Re=De/2;
prompt10='Enter initial parabola angle: ';
prompt11='Enter final parabola angle: ';
Thetaip=input(prompt10);
Thetafp=input(prompt11);
if isempty(Thetaip)
    Thetaip=36
end
if isempty(Thetafp)
    Thetafp=13
end
figure(2)
circled(0, 2.5*Rt, 1.5*Rt, [200:.1:270]);
circled(0, 1.382*Rt, 0.382*Rt, [270:.1:(270+Thetaip)]);
hold on
axis equal
ax=gca;
ax.XTick=[-1 -.5 0 .5 1 1.5 2 2.5];
ax.YTick=[-1 -.75 -.5 -.25 0 .25 .5 .75 1];
xlim([-.5,2.5]);
ylim([-1,1]);
xfit=[sind(36)*0.382*Rt;(sind(36)*0.382*Rt)+.00001;Lb;Lb+.01];
yfit=[Rt+(sind(90)-sind(36))*0.382*Rt*.5;Rt+.00001+(sind(90)-sind(36))*0.382*Rt*.5;Re;Re+.01*sind(Thetafp)];
yay=polyfit(xfit,yfit,3);
yay1=polyval(yay,[sind(36)*0.382*Rt:.01:Lb]);
plot([sind(36)*0.382*Rt:.01:Lb],yay1,'b')
plot([sind(36)*0.382*Rt:.01:Lb],-yay1,'b')
grid on

hold off
f1=-(cwall/2)/Lb;
f2=(Re-Rt)/Lb;
Mne=2*pi*rho_c*Lb*(3^-1*f1*f2*Lb^2+.5*((f1*Rt)+(f2*cwall*.5))*Lb + (Rt*cwall*.5));
Mne=ceil(Mne)         %You know,
Mn=2*Mne             %For safety

Mengt=(Mc+Mn)/.4;
Minj=.249*Mengt;
Mabl=.352*Mengt;
Mtcp=Meng-(Mc+Mn);


%% Tank mass
prompt11='Choose a tank mass factor: ';
PhiTank=input(prompt11);
if isempty(PhiTank)
PhiTank=2500
end
prompt12='Press enter for Helium or type pressurant name: '
pressurant=input(prompt12);
if isempty(pressurant);
    kpres=1.66;
    MMpres=4.003;
else
    prompta="Enter pressurant's gamma value: ";
    promptb="Enter pressurant's molecular mass: ";
    kpres=input(prompta);
    MMpres=input(promptb);
end
M_RP1_tank=2*P_ftank*Vt_f/(PhiTank*g);
M_o_tank=2*P_otank*Vt_ox/(PhiTank*g);
Pave=(P_ftank+P_otank)/2;
P_He=21000000;
Tini=273;
Tfin=Tini*(Pave/P_He)^((kpres-1)/kpres);


 
Mpres=1.05*Pave*(Vt_f+Vt_ox)*MMpres/(8314*Tfin);
Vpres=0;
Vpres2=1;
while (Vpres~=Vpres2)
Vpres=Mpres*8312*Tini/(P_He*MMpres);
Mpres2=1.05*Pave*(Vt_f+Vt_ox+Vpres)*MMpres/(8314*Tfin);
Vpres2=Mpres2*8312*Tini/(P_He*MMpres);
Mpres=1.05*Pave*(Vt_f+Vt_ox+Vpres2)*MMpres/(8314*Tfin);
end
Mprest=P_He*Vpres/(g*6350);
Mtvc=Meng-Mc-Mne;
Mstr=(Mc+Mn+M_o_tank+M_RP1_tank+Mprest)*.1;
Mtot=Mprest+Mpres+M_o_tank+Mo+M_RP1_tank+Mf+Mengt+Mtcp+Mstr;
tburn=(Mo+Mf)/mdot;
Itot=Isp*(Mo+Mf);
Mthr=Mc+Mn;
odot=2.3*(mdot/3.3);
fdot=mdot/3.3;
table(F,tburn,Isp,Itot,mdot,fdot,odot)
table(Mf,Mo,Mpres,Mthr,M_RP1_tank,M_o_tank,Mprest,Mtcp,Mstr)

function h = circled(x,y,r,th)
hold on
xunit = r * cosd(th) + x;
yunit = r * sind(th) + y;
h = plot(xunit, yunit, 'b');
h = plot(xunit,-yunit, 'b');
hold off
end
