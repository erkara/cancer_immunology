% replicate figure5a treatment
% 3 days of chemo + 2 day of rest + cart injection on day 5
% chemo: high strength for one half hour each day for 3 days (C3), 125
% cart: 1e7 at day 5
clc;clear all;close all;

%figure name to save and cart dose
figname = "fig5awtreat.jpg";

%time parameters
dt = 0.01;T_max = 100;tspan = 0:dt:T_max;
%chemo injection days
c1 = 1;c2 = 2;c3 = 3;
index1 = find(tspan==c1);index2 = find(tspan==c2);index3 = find(tspan==c3);
%cart injecttion day(s)
c4 = 5;
index4 = find(tspan==c4);
chemo = 5;
cart_dose = 1e7;
%solve ode system until day1
t1 = tspan(1:index1);
T0 = 1e10;
E0 = 1e5;
CT0 = 0;
M0 = chemo;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t1,[T0,E0,CT0,M0]);
T = y(:,1);
E = y(:,2);
CT = y(:,3);
M = y(:,4);
%solve ode system until from day1 to day2
t2 = tspan(index1+1:index2);
T0 = T(end);
E0 = E(end);
CT0 = CT(end);
M0 = M(end) + chemo;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t2,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];
%solve ode system until from day2 to day3
t3 = tspan(index2+1:index3);
T0 = T(end);
E0 = E(end);
CT0 = CT(end);
M0 = M(end) + chemo;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t3,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];
% we are done with chemo, now 1 day of rest, just one ode
t4 = tspan(index3+1:index4);%days from 3 to 4
T0 = T(end);
E0 = E(end);
CT0 = CT(end);
M0 = M(end);
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t4,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];
%now cart injection
t5 = tspan(index4+1:end);%day5
T0 = T(end);
E0 = E(end);
CT0 = cart_dose*0.01;
M0 = M(end);
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t5,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];

%inspect the results
%plot(tspan,M)
%xlim([0,4])

semilogy(tspan,E,"b-",LineWidth=2)
hold on
semilogy(tspan,T,"r-",LineWidth=2)
hold on
semilogy(tspan,CT,color="#EDB120",LineWidth=2)
legend("Effector Cells","Tumor Cells","Cart Cells")
ylim([1,1e15])
xlabel('Time(days)')
ylabel('Cells')
% saveas(gcf,'fig3b_treatment.png')




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

    dydt(4) = -gamma * M;%M->y(4)
end