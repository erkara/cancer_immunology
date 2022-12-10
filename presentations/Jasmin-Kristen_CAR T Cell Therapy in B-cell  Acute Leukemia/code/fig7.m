%figure7 results 
%FIGURE OUT HOW TO CHANGE SCALE TO MONTHS
clc;clearvars;
dt = 0.1;
 
T_max = 900;
tspan = 0:dt:T_max;
 
C0=10e7;
L0=5e10;
B0=2.5e10;
[~,y] = ode23s(@(t,y) Case1(t,y),tspan,[C0,L0,B0]);
 
%store all variables
CarT = y(:,1);
Leuk = y(:,2);
Bc = y(:,3);
 
plot(tspan,Leuk,'r-',linewidth=2)
hold on 
plot(tspan,CarT,'g-',linewidth=2)
hold on
plot(tspan,Bc,'b-',linewidth=2)
legend({'Leukemia','CAR T','B'})
grid on
ylim([0,8e10])
xlabel('Time (days)') 
ylabel('Cell Number (x10^10)')
title('Figure 7')
figure;

[~,y] = ode23s(@(t,y) Case2(t,y),tspan,[C0,L0,B0]);
 
%store all variables
CarT2 = y(:,1);
Leuk2 = y(:,2);
Bc2 = y(:,3);
 
plot(tspan,Leuk2,'r-',linewidth=2)
hold on
plot(tspan,CarT2,'g-',linewidth=2)
hold on
plot(tspan,Bc2,'b-',linewidth=2)
legend({'Leukemia','CAR T','B'})
grid on
ylim([0,8e10])
xlabel('Time (days)') 
ylabel('Cell Number (x10^10)')
title('Figure 7')
figure;

[~,y] = ode23s(@(t,y) Case3(t,y),tspan,[C0,L0,B0]);
 
%store all variables
CarT3 = y(:,1);
Leuk3 = y(:,2);
Bc3 = y(:,3);
 
plot(tspan,Leuk3,'r-',linewidth=2)
hold on
plot(tspan,CarT3,'g-',linewidth=2)
hold on
plot(tspan,Bc3,'b-',linewidth=2)
legend({'Leukemia','CAR T','B'})
grid on
ylim([0,8e10])
xlabel('Time (days)') 
ylabel('Cell Number (x10^10)')
title('Figure 7')
figure;
 
function dydt = Case1(t,y)
    dydt = zeros(3,1);
    %parameters from the paper
    tau_b = 60;
    rho_l = 1/30;
    tau_c = 14;
    alpha = 4.5e-11;
    rho_c = .25*alpha;
    car50 = 10^9;
    inj0 = 10^7;
    tau_i = 6;
    beta = 0.1;

    dydt(1) = rho_c*(y(2)+y(3))*y(1) + ((rho_c*beta*inj0)/(1+(y(1)/car50)))*y(1) - (y(1)/tau_c);
    dydt(2) = rho_l*y(2) - alpha*y(2)*y(1);
    dydt(3) = ((inj0/tau_i)/(1+(y(1)/car50))) - alpha*y(3)*y(1) - (y(3)/tau_b);
end
function dydt = Case2(t,y)
    dydt = zeros(3,1);
    %parameters from the paper
    tau_b = 60;
    rho_l = 1/30;
    tau_c = 14;
    alpha = 4.5e-11;
    rho_c = .25*alpha;
    car50 = 10^9;
    inj0 = 3e8;
    tau_i = 6;
    beta = 0.1;

    dydt(1) = rho_c*(y(2)+y(3))*y(1) + ((rho_c*beta*inj0)/(1+(y(1)/car50)))*y(1) - (y(1)/tau_c);
    dydt(2) = rho_l*y(2) - alpha*y(2)*y(1);
    dydt(3) = ((inj0/tau_i)/(1+(y(1)/car50))) - alpha*y(3)*y(1) - (y(3)/tau_b);
end
function dydt = Case3(t,y)
    dydt = zeros(3,1);
    %parameters from the paper
    tau_b = 60;
    rho_l = 1/30;
    tau_c = 14;
    alpha = 4.5e-11;
    rho_c = .25*alpha;
    car50 = 10^9;
    inj0 = 5e8;
    tau_i = 6;
    beta = 0.1;

    dydt(1) = rho_c*(y(2)+y(3))*y(1) + ((rho_c*beta*inj0)/(1+(y(1)/car50)))*y(1) - (y(1)/tau_c);
    dydt(2) = rho_l*y(2) - alpha*y(2)*y(1);
    dydt(3) = ((inj0/tau_i)/(1+(y(1)/car50))) - alpha*y(3)*y(1) - (y(3)/tau_b);
end