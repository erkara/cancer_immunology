%figure6b: due to randomness pointed out in virtual patient
%generation process, we are not getting the same results.
clc;clearvars;close all;
format bank;
data = table2array(readtable('results/vp.csv'));
s = data(:,1:9);
%take a subset
s = s(1:100,:);
[n,l] = size(s);
%doses = [3e6,2.5e6,2e6,1.5e6,1e6,0.5e6];
doses = [3e6,2.5e6,1.5e6,1e6,0.5e6,0.1e6];
m = length(doses);
%time parameters
dt = 0.1;T_max = 300;tspan = 0:dt:T_max;
T_inject = 42;
index = find(tspan == T_inject);
t1 = tspan(1:index);
t2 = tspan(index+1:end);
surv_days = zeros(n,m);
progress = waitbar(0, 'Starting');
k = 0;
%loop over patients
for i=1:n
    %loop over doses
    for j=1:m
        k = k+1;
        waitbar(k/(n*m),progress,sprintf('remain = %0.1f %%',100*(1-k/(n*m))));
        T = zeros(length(tspan),1);
        CT0 = 0;
        CM0 = 0;
        T0 = 2e6;
        [~,y] = ode45(@(t,y) GetHDLM(t,y,s(i,:)),t1,[CT0,CM0,T0]);
        T(1:index) = y(:,3);
        
        %day42, start terapy
        CT0 = doses(j);
        CM0 = 0;
        T0 = T(index);

        [~,y] = ode23s(@(t,y) GetHDLM(t,y,s(i,:)),t2,[CT0,CM0,T0]);
        T(index+1:end) = y(:,3);

        tumor_max = T(end);
        thresold = 1e10;
        fprintf("tumor_start: %0.1d tumor_max %0.1d\n",T(index),tumor_max)
        if tumor_max > thresold
            %get the first day exceeding thresold
            surv_days(i,j) = tspan(find(T>thresold,1));
        else
            %full survival
            surv_days(i,j) = 300 ;
        end
    end
    %writematrix([s(i,:),surv_days(i,:)],'results/doses_vs_surdays.csv','WriteMode','append');
end
delete(progress)
%plot the results
d = 1:300;
for i=1:length(doses)
    percents = zeros(length(d),1);
        for j=1:length(d)
            percents(j) = round(100*(length(find(surv_days(:,i) > d(j))))/n,2); 
        end
    survival = round((length(find(surv_days(:,i)==300))/n)*100,1);
    plot(d,percents,'DisplayName',sprintf('dose: %0.1e\nsurvial: %0.2f %%',...
        doses(i),survival),linewidth = 1.5);
    grid on;
    legend("Location","southwest");
    hold on;
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


