%Generate figure4-def results, fractional dose strategy
%fig4d: keep the default doses and days
%fig4e: doses = [0.5e6,0.5e6,0.5e6,0.5e6];days = [42,56,70,84]; 
%fig4f: doses = [0.1,0.3,0.6]*2e6,days = [42,43,44];
clc;close all;clearvars;

%time parameters
dt = 0.01;Tmax = 200;tspan = 0:dt:Tmax;
%modify here based on the dosing schedule.
doses = [0.5e6,0.5e6,0.5e6,0.5e6];
days = [42,49,56,63];
%cut-off indices for tspan.we will solve our ODE between these times
c = ismember(tspan,days);
indices = find(c);
indices=[1,indices,length(tspan)];
%tumor grows until day-42
CT0 = 0;
CM0 = 0;
T0 = 2e6;
k = 0;
temp_indices = indices(1)+k:indices(2);
[~,y] = ode23s(@(t,y) GetHDLM(t,y),tspan(temp_indices),[CT0,CM0,T0]);
CT(temp_indices) = y(:,1);
CM(temp_indices) = y(:,2);
T(temp_indices) = y(:,3);

for i=1:length(doses)
    k = 1;
    temp_indices = indices(i+1)+k:indices(i+2);
    CT0 = CT(end)+doses(i);
    CM0 = CM(end);
    T0 = T(end);
    [~,y] = ode23s(@(t,y) GetHDLM(t,y),tspan(temp_indices),[CT0,CM0,T0]);
    CT(temp_indices) = y(:,1);
    CM(temp_indices) = y(:,2);
    T(temp_indices) = y(:,3);
end


yyaxis left
plot(tspan,CT,'g-',linewidth=1)
hold on
plot(tspan,CM,'b-',linewidth=1)
ylim([0,3e6]);
hold on

yyaxis right
plot(tspan,T,'r-',linewidth=2)
legend({'C_T: Effector CART','C_M: Memory CART','T: Tumor'},'Location','northwest')
ylim([0,2.5e7])
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
    
    dydt(1) = phi*y(1) - rho*y(1) + theta*y(3)*y(2) - alpha*y(3)*y(1);
    dydt(2) = eps*y(1) - theta*y(3)*y(2) - mu*y(2);
    dydt(3) = r*y(3)*(1 - b*y(3)) - gamma*y(1)*y(3);
end
