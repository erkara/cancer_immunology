%use this file to find the best parameters fitting to an ode or ode-system
%based on existing data. Notice how we build our ode with a free parameter
% "s" to be used in optimization process. You can easily modify this code
%to deploy it for your own problem.
clearvars;clc;close all;
%time-parameters
dt = 0.1;time_max = 30;tspan = 3:dt:time_max;

%a better practice would be to read the data from a file.
days = [3,8,11,15,18,23];
volume =[50,83.4,155.2,271.6,400.3,880.5]';

%lower-upper bound for r and K in order
lb = [0,1e3];
ub = [1,1e4];
init_val = [0.1,2e3];
[optim_sol,sumsq] = FitParam(lb,ub,init_val,tspan,days,volume);

%inspect the results
fprintf('r=%0.3f K=%0.2f ',optim_sol.s(1),optim_sol.s(2));
fprintf(' MAE: %0.2f\n',sumsq);

%plot the optimized model
[tumor_total,Ts] = SolveODE(optim_sol.s,tspan,days);
plot(days,tumor_total,'r-')
hold on;
%plot the real data
plot(days,volume,'bo','MarkerSize',12)
legend({'model','raw-data'},"Location","northwest")
grid on;


function [optim_sol,sumsq] = FitParam(lb,ub,init_val,tspan,days,volume)
    %create an optim variable with certain upper and lower bounds
    s = optimvar('s',length(ub),"LowerBound",lb,"UpperBound",ub);
    %initial values for the variables to be optimized
    s0.s = init_val;
    %solve ODE with optim parameters
    [tumor_total,~] = fcn2optimexpr(@SolveODE,s,tspan,days);
    %optimize this
    obj = sum(sqrt((tumor_total - volume).^2))/length(volume);
    prob = optimproblem("Objective",obj);
    options = optimoptions(prob);
    options.ConstraintTolerance = 1e-15;
    options.Display = 'iter';
    %actual optimization stage. optim_sol holds the optimized variables
    %while sumsq is the value of obj we specified. 
    [optim_sol,sumsq] = solve(prob,s0,'Options',options);
end



function [tumor_total,Ts] = SolveODE(s,tspan,days)
    y0 = 50;
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
    %find indices of the corresponding days
    [~,indices] = ismember(days,tspan);
    %solve ODE with parameter s
    [~,y] = ode45(@(t,y) GetODE(t,y,s),tspan,y0,opts);
    %tumor volume
    Ts = y(:,1);
    %get the tumor volume at the correspoding dates
    tumor_total = Ts(indices);
end
 

%logistic growth model
function dydt = GetODE(t,y,s)
    r = s(1);
    K = s(2);
    dydt = r*y.*(1-y/K);
end
