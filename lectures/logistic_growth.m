%logistic growth model
clc;clearvars;close all;
dt = 0.1;T = 30;tspan = 0:dt:T;
y0 = 2;
r = 0.2;
K = 10;
%solve IVP
[t,y] = ode45(@(t,y) GetODE(t,y,r,K),tspan,y0);
plot(t,y,'r')
%let's define and plot the exact solution
legend('numerical solution')
ylim([0,K+5])
grid on;

function dydt = GetODE(t,y,r,K)
    %this tells us we have a single ODE
    dydt = zeros(1);
    %notice that we use y(1) instead of y since the solution
    %is 1D vector
    dydt(1) = r*y(1)*(1 - y(1)/K);
end