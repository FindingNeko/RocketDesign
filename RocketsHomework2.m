%Rocket Propulsion Homework 2 problem 6
%Bryan Rathke

%Given

Ar=2.3;
Ari=3.02;
At=5;
k=1.3;
R=66;
p3=300;
p2=200;
p1=100;
T0=5300;
pa=10;
g0=32.174;
kp1=k+1;
km1=k-1;

%Equations

    %Pressure ratio between chamber and atmosphere
    r3=p3/pa;
    r2=p2/pa;
    r1=p1/pa;
    %Effective exhaust velocity for given area ratio
    c1=sqrt(((2*k)/km1)*R*T0*(1-(p1/(12*p1))^(km1/k)));
    c2=sqrt(((2*k)/km1)*R*T0*(1-(p2/(12*p2))^(km1/k)));
    c3=sqrt(((2*k)/km1)*R*T0*(1-(p3/(12*p3))^(km1/k)));
    %Ideal exhaust velocity for optimum and actual area ratio
    
    %Propellant flow
    mdot3=At*p3*k*sqrt((2/kp1)^(kp1/km1)/k*R*T0);
    mdot2=At*p2*k*sqrt((2/kp1)^(kp1/km1)/k*R*T0);
    mdot1=At*p1*k*sqrt((2/kp1)^(kp1/km1)/k*R*T0);
    %Thrust
    F1=At*p1*sqrt(2*k^2*(2/kp1)^(kp1/km1)*(1-(pa/p1)^(km1/k)/km1));
    F1=At*p2*sqrt(2*k^2*(2/kp1)^(kp1/km1)*(1-(pa/p2)^(km1/k)/km1));
    F1=At*p3*sqrt(2*k^2*(2/kp1)^(kp1/km1)*(1-(pa/p3)^(km1/k)/km1));
    %Specific impulse
    Is1=F1/(mdot1*g0);
    Is2=F2/(mdot2*g0);
    Is1=F3/(mdot3*g0);
    %Exit pressure
    
    
    %Plots
    
    figure(1);
    subplot(3,5,6)
    plot(p1,pa)
    plot(p2,pa)
    plot(p3,pa)
    