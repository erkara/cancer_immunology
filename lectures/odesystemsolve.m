%generic ode systems solver
clc;clearvars;close all
dt = 0.01;T = 1;tspan = 0:dt:T;
%x1(0)=2,x2(0)=4, the order is important
y0 = [2,4];
%solve IVP
[t,y] = ode45(@(t,y) GetODE(t,y),tspan,y0);
plot(t,y(:,1),'r-')
hold on
plot(t,y(:,2),'b--')
legend('x1(t)','x2(t)')
grid on;
%plot x1 vs x2 on the same plane 
figure;
plot(y(:,1),y(:,2))
grid on;

function dydt = GetODE(t,y)
    %this tells us we have two ODES
    dydt = zeros(2,1);
    %notice that we use y(1),y(2) instead of y since the solution
    %is 2D vector
    dydt(1) = y(1) + 3*y(2) + 12*t - 11;
    dydt(2) = 5*y(1) + 3*y(2) - 3;
end