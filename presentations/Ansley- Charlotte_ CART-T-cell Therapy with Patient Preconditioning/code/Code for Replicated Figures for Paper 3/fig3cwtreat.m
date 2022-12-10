% replicate figure3c treatment
% 5 days of chemo + 2 days rest + cart injection on day 8
% chemo: medium strength for one half hour each day for 5 days (C5), 75
% cart: 10e7 at day 8
clc;clear all;close all;

%figure name to save and cart dose
figname = "fig3ctreat.jpg";

%time parameters
dt = 0.01;T_max = 100;tspan = 0:dt:T_max;
%chemo injection days
c1 = 1;c2 = 2;c3 = 3;c4= 4;c5 = 5;
index1 = find(tspan==c1);index2 = find(tspan==c2);index3 = find(tspan==c3);
index4 = find(tspan==c4);index5 = find(tspan==c5);
%cart injection day(s)
c6 = 8;
index6 = find(tspan==c6);
chemo = 5;
cart_dose = 1e7;
%solve ode system until day1
t1 = tspan(1:index1);
T0 = 1e10;
E0 = 4e5;
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
%solve ode system until from day3 to day4
t4 = tspan(index3+1:index4);
T0 = T(end);
E0 = E(end);
CT0 = CT(end);
M0 = M(end) + chemo;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t4,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];
%solve ode system until from day4 to day5
t5 = tspan(index4+1:index5);
T0 = T(end);
E0 = E(end);
CT0 = CT(end);
M0 = M(end) + chemo;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t5,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];
% we are done with chemo, now 2 days rest, just one ode
t6 = tspan(index5+1:index6);%days from 6 to 7
T0 = T(end);
E0 = E(end);
CT0 = CT(end);
M0 = M(end);
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t6,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];
%now cart injecttion
t7 = tspan(index6+1:end);%day8
T0 = T(end);
E0 = E(end);
CT0 = cart_dose*0.01;
M0 = M(end);
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t7,[T0,E0,CT0,M0]);
T = [T;y(:,1)];
E = [E;y(:,2)];
CT = [CT;y(:,3)];
M = [M;y(:,4)];

%inspect the results
%plot(tspan,M)
%xlim([0,6])

semilogy(tspan,E,"b-",LineWidth=2)
hold on
semilogy(tspan,T,"r-",LineWidth=2)
hold on
semilogy(tspan,CT,color="#EDB120",LineWidth=2)
legend("Effector Cells","Tumor Cells","Cart Cells")
ylim([1,1e15])
% saveas(gcf,'fig3b_treatment.png')




function dydt = GetHDLM(~,y)
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

    dydt(4) = -gamma * M;%M->y(4)
end