%solve lotka-volterra model
clc;clearvars;close all;
%time parameters
dt = 0.1;T = 1000;tspan = 0:dt:T;
%x1(0)=80,x2(0)=50, the order is important
y0 = [80,50];
%solve IVP
a = 0.2;
b = 0.0025;
c = 0.01;
d = 0.002;
[t,y] = ode45(@(t,y) GetODE(t,y,a,b,c,d),tspan,y0);
plot(t,y(:,1),'r-')
hold on
plot(t,y(:,2),'b-')
legend('prey','predator','FontSize',10)
grid on;
%plot x1-->x axis vs x2--yaxis 
figure;
plot(y(:,1),y(:,2))
xlabel('prey','FontSize',10)
ylabel('predator','FontSize',10)
grid on;

function dydt = GetODE(t,y,a,b,c,d)
    %this tells us we have two ODES
    dydt = zeros(2,1);
    %notice that we use x1-->y(1),x2-->y(2) 
    dydt(1) = a*y(1) - b*y(1)*y(2);
    dydt(2) = -c*y(2) + d*y(1)*y(2);
end