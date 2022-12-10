%figure3a results 
clc;clearvars;
dt = 0.1;
 
T_max = 60;
tspan = 0:dt:T_max;

%initial conditions from figure 2
C0=5e5;
L0=5e10;
B0=2.5e10;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),tspan,[C0,L0,B0]);
 
%store all solid line variables
CarT = y(:,1);
Leuk = y(:,2);
Bc = y(:,3);

C0=15e5;
L0=5e10;
B0=2.5e10;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),tspan,[C0,L0,B0]);
 
%store all dashed variables
CarTD = y(:,1);
LeukD = y(:,2);
BcD = y(:,3);

C0=45e5;
L0=5e10;
B0=2.5e10;
[~,y] = ode23s(@(t,y) GetHDLM(t,y),tspan,[C0,L0,B0]);
 
%store all dotted variables
CarTDo = y(:,1);
LeukDo = y(:,2);
BcDo = y(:,3);
 
plot(tspan,Leuk,'r-',linewidth=2)
hold on 
plot(tspan,CarT,'g-',linewidth=2)
hold on
plot(tspan,Bc,'b-',linewidth=2)
hold on
plot(tspan,LeukD,'r--',linewidth=2)
hold on 
plot(tspan,CarTD,'g--',linewidth=2)
hold on
plot(tspan,BcD,'b--',linewidth=2)
hold on
plot(tspan,LeukDo,'r:',linewidth=2)
hold on 
plot(tspan,CarTDo,'g:',linewidth=2)
hold on
plot(tspan,BcDo,'b:',linewidth=2)
legend({'Leukemia','CAR T','B'})
grid on
ylim([0,6e10])
xlabel('Time (days)') 
ylabel('Cell Number (x10^10)')
title('Figure 3a')
 
function dydt = GetHDLM(t,y)
    dydt = zeros(3,1);
%parameters from the paper
    tau_b = 60;
    rho_l = 1/30;
    tau_c = 14;
    alpha = 4.5e-11;
    rho_c = .25*alpha;

    dydt(1) = rho_c*(y(2)+y(3))*y(1) -(1./tau_c)*y(1);
    dydt(2) = rho_l*y(2) - alpha*y(2)*y(1);
    dydt(3) = -alpha*y(3)*y(1) - (1/tau_b)*y(3);
end