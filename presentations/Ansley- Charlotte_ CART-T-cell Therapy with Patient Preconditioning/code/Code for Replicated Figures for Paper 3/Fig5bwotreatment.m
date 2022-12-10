% replicate figure5b no treatment
clc;clear all;close all;

%figure name to save and cart dose
figname = "fig5bnotreat.jpg";


%time parameters
dt = 0.01;T_max = 100;tspan = 0:dt:T_max;

T0 = 5e9;
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
legend("Effector Cells no treatment","Tumor Cells no treatment")
ylim([1e0,1e15])
xlabel('Time(days)')
ylabel('Cells')
%saveas(gcf,'fig3bnotreat.png')



function dydt = GetHDLM(t,y)
    dydt = zeros(4,1);
    %parameters from the paper(patient 1)
    a = 2.55e-1;
    b = 2e-12;
    dE = 2.03;
    dC = 2.25;
    g = 1.4e3;
    jE = 1.1e-2;
    jC = 2.42e-1;
    K = 1.65e9;
    k = 2.019e5;
    l = 1.395;
    mE = 7e-3;
    mC = 2.93e-2;
    qE = 3.42e-11;
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

