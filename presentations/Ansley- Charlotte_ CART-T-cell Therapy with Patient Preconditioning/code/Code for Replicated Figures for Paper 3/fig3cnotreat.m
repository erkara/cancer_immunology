% replicate figure3c no treatment
clc;clear all;close all;

%figure name to save and cart dose
figname = "fig3cnotreat.jpg";


%time parameters
dt = 0.01;T_max = 50;tspan = 0:dt:T_max;

T0 = 1e10;
E0 = 1e5;
CT0 = 0;
M0 = 0;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),tspan,[T0,E0,CT0,M0]);

T = y(:,1);
E = y(:,2);
CT = y(:,3);
M = y(:,4);

%plot the results.  

semilogy(tspan,E,"b-")
hold on
plot(tspan,T,"r--")
legend("Effector Cells","Tumor Cells")
ylim([1e0,1e15])
saveas(gcf,'fig3bnotreat.png')



function dydt = GetHDLM(t,y)
    dydt = zeros(4,1);
    %parameters from the paper(patient 2)
    a = 1.76e-1;
    b = 2e-12;
    dE = 2.03;
    dC = 2.25;
    g = 4.7e4;
    jE = 7.46e-3;
    jC = 1.65e-1;
    K = 1.65e9;
    k = 7.0e7;
    l = 1.419;
    mE = 3.4e-2;
    mC = 2.93e-2;
    qE =6.71e-11;
    qC = 3.0e-11;
    s = 3.05e-1;
    KT = 7.00e-1;
    KE = 6.00e-1;
    KC = 6.00e-1;
    gamma = 9.00e-1;
    T = y(1);
    E = y(2);
    C = y(3);
    M = y(4);

    %ODE systems
    r1 = (E/T)^l;
    DE = dE*(r1/(s+r1))*T;
    r2 = (C/T)^l;
    DC = dC*(r2/(s+r2))*T;

    dydt(1) = a*T*(1 - b*T) - DE - DC - KT*(1 - exp(-M))*T;%T->y(1)

    dydt(2) = g - mE*E - jE*log((E+C)/K)*((DE^2)/(k+DE^2))*E - qE*E*T - KE*(1-exp(-M))*E;%E->y(2)

    dydt(3) = -mC*C-jC*log((E+C)/K)*((DC^2)/(k+DC^2))*C - qC*C*T - KC*(1-exp(-M))*C;%C->y(3)

    dydt(4) = - gamma * M;%M->y(4)
end
