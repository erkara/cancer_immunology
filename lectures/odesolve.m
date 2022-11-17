dt = 0.1;
T = 5;
tspan = 0:dt:T;
y0 = 4;
[t,y] = ode45(@(t,y) GetODE(t,y),tspan,y0);
plot(t,y,'r')
%let's define and plot the exact solution
exact = @(t) t.^2+4;
hold on
plot(t,exact(t),'b--');
legend('numerical solution','exact solution')




function dydt = GetODE(t,y)
    %this tells us we have a single ODE
    dydt = zeros(1);
    %notice that we use y(1) instead of y since the solution
    %is 1D vector
    dydt(1) = 2*t*y(1)/(t^2+4);
end