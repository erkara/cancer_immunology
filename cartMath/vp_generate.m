%use this code to generate virtual patients as described in the paper.
%carefully examine our acceptance criteria in "if" statement below.
%However, we diverge from the paper at this point since this process
%involves some randomness. Run fig6a.m and observe that percetages
%are little bit off comparing to the paper.
clc;clearvars;close all;
dt = 0.1;
k = 0;
T_inject = 42;
T_max = 300;
tspan = 0:dt:T_max;
index = find(tspan == T_inject);
t1 = tspan(1:index);
t2 = tspan(index+1:end);
N = 80000 ;
VP = zeros(N,1);
min_sample = 4000;
progress = waitbar(0,'starting...');
tic
for i = 1:N
    waitbar(i/N,progress,sprintf('%0.2f %%',100*(1-i/N)));
    T = zeros(length(tspan),1);
    s = GetUniformSample();
    %tumor grows without intervention
    C_T0 = 0.;
    C_M0 = 0.;
    T0 = 2e6;
    [~,y] = ode45(@(t,y) GetHDLM(t,y,s),t1,[C_T0,C_M0,T0]);
    T(1:index) = y(:,3);
    %t = 42 days, inject 2e6 CART-Cells  
    C_T0 = 2e6;
    C_M0 = 0;
    T0 = max(T);
    [~,y] = ode23s(@(t,y) GetHDLM(t,y,s),t2,[C_T0,C_M0,T0]);
    %get the max tumor burden
    T(index+1:index+length(y(:,3))) = y(:,3);
    tumor_max = max(T);
    thresold = 1e10;
   
    if tumor_max > thresold
        fprintf("r = %0.3f\n",s(7))
        day_index = find(T>thresold,1);
        k = k + 1;
        if mod(k,10)==0
            fprintf('k=%i i=%i occurance=%0.3f \n',k,i,k/i)
            toc
        end
        %the day thresold first occured and parameters for this case 
        VP(k) = tspan(day_index);
        %writematrix(tspan(day_index),'thresold_days.csv','WriteMode','append');
        writematrix([s,tspan(day_index)],'results/vp.csv','WriteMode','append');
        %median survival of all patient VPs so far
        med_survival = median(nonzeros(VP));
        %we want to collect more than "min_sample"
        if k>min_sample
            if med_survival<128+10 && med_survival>128-10
                fprintf('desired constraint found at %i',k)
                    break
            end
        end
    end
end



 
function dydt = GetHDLM(t,y,s)
    %[phi, rho, eps,  theta, alpha, mu,  r,    b,  gamma] = s;
    phi = s(1); 
    rho = s(2);
    eps = s(3);  
    theta = s(4); 
    alpha = s(5);
    mu = s(6);
    r =  s(7);    
    b = s(8);
    gamma = s(9);
    %s(1),s(2),s(3), s(4),  s(5), s(6), s(7),s(8),s(9)
    dydt = zeros(3,1);
    dydt(1) = phi*y(1) - rho*y(1) + theta*y(3)*y(2) - alpha*y(3)*y(1);
    dydt(2) = eps*y(1) - theta*y(3)*y(2) - mu*y(2);
    dydt(3) = r*y(3)*(1 - b*y(3)) - gamma*y(1)*y(3);
end


function s = GetUniformSample()
    %s =  [phi,rho,eps,theta,alpha,mu,r,b,gamma]
    phi0 = 0.265; 
    rho0 = 0.350;
    eps0 = 0.150;
    theta0 = 6.0e-6;
    alpha0 = 4.5e-8;
    mu0 = 5.0e-3;
    r0 = 5.650026e-2;
    b0 = 1.404029e-12;
    gamma0 = 3.715843e-6;
    param_array = [phi0,rho0,eps0,theta0,alpha0,mu0,r0,b0,gamma0];
    a = param_array - 0.6*param_array;
    b = param_array + 0.6*param_array;
    s = a+(b-a).*rand(length(param_array),1)';
    %phi<rho<phi+eps  
    s(2) = s(1)+ s(3)*rand();
end