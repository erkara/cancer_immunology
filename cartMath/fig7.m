%figure7 data: this code generates the data required for sensitivity analysis
%described in figure7. We could not generate the same results with
% 1e6 CART cells but we are getting similar results for CR cases with 2.5e6 
% as we discussed in the class. To produce the actual figures, you need to run the Jupyter
%notebook "sensitivity.ipynb". 

clc;clearvars;close all;format bank;
data = table2array(readtable('results/vp.csv'));
s = data(:,1:9);
%s = s(1:100,:);
[N,l] = size(s);
dt = 0.1;T_max = 300;tspan = 0:dt:T_max;
T_inject = 42;
index = find(tspan == T_inject);
t1 = tspan(1:index);
t2 = tspan(index+1:end);
index55 = find(tspan==55);
index75 = find(tspan==75);
progress = waitbar(0, 'Starting');
k = 0; 
col_names = ["phi","rho","eps","theta","alpha",...
    "mu","r","b","gamma","CM55","CM75","Response"];
for i=1:N
    CT = zeros(length(tspan),1);
    CM = zeros(length(tspan),1);
    T = zeros(length(tspan),1);

    k = k+1;
    waitbar(k/(N),progress,sprintf('remain = %0.1f %%',100*(1-k/(N))));
    %let the tumor grow until day=42
    C_T0 = 0;
    C_M0 = 0;
    T0 = 2e6;
    [~,y] = ode45(@(t,y) GetHDLM(t,y,s(i,:)),t1,[C_T0,C_M0,T0]);
    CT(1:index) = y(:,1);
    CM(1:index) = y(:,2);
    T(1:index) = y(:,3);
    
    %day42, start terapy
    C_T0 = 2.5e6;
    C_M0 = 0;
    T0 = max(T);
    a = T0;
    

    [~,y] = ode23s(@(t,y) GetHDLM(t,y,s(i,:)),t2,[C_T0,C_M0,T0]);
    CT(index+1:end) = y(:,1);
    CM(index+1:end) = y(:,2);
    T(index+1:end) = y(:,3);
    


    %get reference values
    CM55 = CM(index55);
    CM75 = CM(index75);
    thresold = 1e10;
    
    if max(T)>thresold
        c = "NR";
    else
        c = "CR";
    end
    data = [s(i,:),CM55,CM75,c];
    data = array2table(data,"VariableNames",col_names);
    writetable(data,'results/fig7_data.csv',...
        'WriteMode','append');
   
end
    
close(progress)






function dydt = GetHDLM(t,y,s)
    dydt = zeros(3,1);
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
    dydt(1) = phi*y(1) - rho*y(1) + theta*y(3)*y(2) - alpha*y(3)*y(1);%CT
    dydt(2) = eps*y(1) - theta*y(3)*y(2) - mu*y(2);%CM
    dydt(3) = r*y(3)*(1 - b*y(3)) - gamma*y(1)*y(3);%T
end


