%HDLM-2 cell-line with CAR-T 123 therapy.
%figure3a results
clc;clearvars;close all;
%time paramaters
dt = 0.1;T_max = 500;tspan = 0:dt:T_max;
T_inject1 = 42;
T_inject2 = 250;
index1 = find(tspan == T_inject1);
index2 = find(tspan == T_inject2);
%tumor grows without intervention
t1 = tspan(1:index1);
CT0 = 0.;
CM0 = 0.;
T0 = 2e6;
[~,y] = ode45(@(t,y) GetHDLM(t,y),t1,[CT0,CM0,T0]);
CT = y(:,1);
CM = y(:,2);
T = y(:,3);
%t = 42 days, inject 2e6 CART-Cells  
t2 = tspan(index1+1:index2);
CT0 = 2e6;
CM0 = CM(end);
T0 = T(end);
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t2,[CT0,CM0,T0]);
CT = [CT;y(:,1)];
CM = [CM;y(:,2)];
T = [T;y(:,3)];

%t = 250 days, challange with 1e6 tumor cells  
t3 = tspan(index2+1:end);
CT0 = CT(end);
CM0 = CM(end);
T0 = T(end) + 1e6;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),t3,[CT0,CM0,T0]);
CT = [CT;y(:,1)];
CM = [CM; y(:,2)];
T = [T;y(:,3)];

%plot the results. Zoom on T to see how tumor cells are eliminated quickly 
yyaxis left
plot(tspan,CT,'g-',linewidth=1)
hold on
plot(tspan,CM,'b-',linewidth=1)
hold on

yyaxis right
plot(tspan,T,'r-',linewidth=2)
legend({'C_T: Effector CART','C_M: Memory CART','T: Tumor'},'Location','northwest')
grid on


function dydt = GetHDLM(t,y)
    dydt = zeros(3,1);
    %parameters from the paper
    phi = 0.265;
    rho = 0.350;
    eps = 0.150;
    theta = 6.0e-6;
    alpha = 4.5e-8;
    mu = 5.0e-3;
    r = 5.650026e-2;
    b = 1.404029e-12;
    gamma = 3.715843e-6;
    %ODE systems
    dydt(1) = phi*y(1) - rho*y(1) + theta*y(3)*y(2) - alpha*y(3)*y(1);%CT
    dydt(2) = eps*y(1) - theta*y(3)*y(2) - mu*y(2);%CM
    dydt(3) = r*y(3)*(1 - b*y(3)) - gamma*y(1)*y(3);%T
end
