clc; close all;

data_a = table2array(readtable('Swanson_Fig4a.csv'));
err_data_a = table2array(readtable('Swanson_4a_Final_errors.csv'));
y_err_u = err_data_a(:,2);
y_err_l = err_data_a(:,3);
time_vals_a = data_a(:,1);
tumor_vals_a = data_a(:,2);
data_b = table2array(readtable('Swanson_Fig4b.csv'));
time_vals_b = data_b(:,1);
tumor_vals_b = data_b(:,2);
dt = 0.1;  
T = 60; 
tspan = 0:dt:T;
T0 = 3.4e7;
S0= 0.01*T0;
E0 = 0;
C0=0;
B0=5;
R0=2.5;
 
[~,y] = ode23s(@(t,y) GetODE(t,y),tspan,[S0, T0, E0,C0, B0, R0]);
S = y(:,1);
T = y(:,2);
E = y(:,3);
B = y(:,4);
C = y(:,5);
R = y(:,6);

%Figure 4b
figure;
plot(tspan, S, 'Linewidth', 1.0);
hold on
plot(time_vals_b, tumor_vals_b, 'ks');
ylim([0,4e6]);
xlim([0,60]);
ylabel("Cancer Stem Cells")
xlabel("t(days)")
title('Figure 4b - Swanson 2021')

%Figure 4a
G = T + S;
figure;
plot(tspan, G, 'Linewidth', 1.0);
hold on
errorbar(time_vals_a,tumor_vals_a,y_err_l, y_err_u, 'ks');
ylim([0,8e8]);
xlim([0,60]);
ylabel("Tumor Cells")
xlabel("t(days)")
title('Figure 4a - Swanson 2021')


%function that is called to generate each model value for the ODE

function dydt = GetODE(t,y)
    dydt = zeros(6,1);
    f = 0.0011;
    v0 = 0.3628;
    p = 0.1351;
    sE = 1e-5;
    r = 0.0125;
    b = 0.4843;
    c_th = 0.1224; 
    n = 0.0454;
    uCS = 1e-6;
    uET = 2e-5;
    c_o = 100;
    K = 1e9;
    uRE = 1e-5;
    sR = 1e-5;
    a = 0.5917;
    c_t = 1.875e9;
    alpha = 1e-8; 
    sB = 7e-4;
    uRC = 1e-5;
    sCT = 1e-4;

    S = y(1);
    T = y(2);
    E = y(3);
    C = y(4);
    B = y(5);
    R = y(6);

    dydt(1) = p* v0*(1 - (T+S)/K)*S - uCS*C*S;
    dydt(2) = n * T *(1 - (T+S)/K) + (1-p)*v0*(1 - ((T+S)/K))*S - uET*(T*E/(1+c_o*B));
    dydt(3) = b*E + (f*E*T)/(1+ (c_th *T*B)) - uRE *R*E - sE*E - r*E; 
    dydt(4) = b*C + ((f*C*S)/(1 +(c_th*T*B))) - uRC*R*C - sCT*C;
    dydt(5) = (a* (T^2))/(c_t + (T^2)) + alpha*R - sB*B;
    dydt(6) = r*E - sR*R;
end
