%Figure-5, diagram of occurrence of complete response (CR: T(300) ≤ 8e5 cells, 
% green dots) and non-response (NR: T (300) ≥ 10e10 cells, red dots) 
% for the HDLM-2 + CAR-T 123 scenario
clc;clearvars;close all;
%time params
dt = 0.1;day_max = 300;tspan = 0:dt:day_max;

%initial tumor numbers changing from 100K to 10M
dT = 1e5;     %100K
T_Max = 1e7;  %10M 
T0_vals = 8e5:dT:T_Max;

%CART cells injection numbers changing from 20K TO 2M
dCT = 2e4;
CT_Max = 2e6;
CT0_vals = 1e5:dCT:CT_Max;

%understand this piece
day_inject = 42;
index = find(tspan == day_inject);
t1 = tspan(1:index);
t2 = tspan(index+1:end);
m = length(T0_vals);
n = length(CT0_vals);
size = m*n;
responses = repmat(string(0),size,1);
CART_Dose =zeros(size,1);
Tumor_Init = zeros(size,1);
k = 0;
%this simulation takes quite a bit time, go grab coffe as progress bar
%displays the progress of the simuation
progress = waitbar(0,'starting...');
counter = 0;
for T0 = T0_vals
    %let the tumor grow for 42 days with initial T0 amount
    CT0 = 0;
    CM0 = 0;
    [~,y] = ode15s(@(t,y) GetHDLM(t,y),t1,[CT0,CM0,T0]);
    T_temp = y(end);
    %for a fixed tumor number T_temp at day=42, use different CART doses
    for CT0 = CT0_vals
        counter = counter + 1;
        waitbar(counter/size,progress,sprintf('%0.1f %%',100*(1-counter/size)));
        [~,z] = ode15s(@(t,z) GetHDLM(t,z),t2,[CT0,CM0,T_temp]);
        T_final = z(end);
        k = k+1;
        Tumor_Init(k) = T0;
        CART_Dose(k) = CT0;
        %classify the outcome of treatment based on the criteria above
        if T_final < 8e5
            responses(k) = "CR"; 
        elseif T_final > 1e10
            responses(k) = "NR";
        else
            responses(k) = "PR";
        end
    end
end
delete(progress)
%write to csv just in case
data = table(CART_Dose,Tumor_Init,responses);
writetable(data,'results/dose_response.csv')
%read the data
data = readtable('results/dose_response.csv');
color = 'gr';
markers = 'ooo';
sz = [4,5,4];

%do scatter plot
%gscatter(x,y,g,clr,sym,siz,doleg,xnam,ynam)
gsh = gscatter(data.CART_Dose,data.Tumor_Init,data.responses,color,markers, ...
    sz,'on','C_T(#cells)','T(#cells)');
xlim([6e4,2.1e6])
ylim([6e5,1.01e7])
grid on
for g = 1:length(gsh)
    gsh(g).MarkerFaceColor = gsh(g).Color;
end



function dydt = GetHDLM(t,y)
    dydt = zeros(3,1);
    %parameters from the paper
    phi = 0.265;
    rho = 0.350;
    eps = 0.150;
    theta = 6.0e-6;
    alpha = 4.5e-8;
    mu = 5.0e-3;
    r = 5.650026e-2;
    b = 1.404029e-12;
    gamma = 3.715843e-6;
    
    dydt(1) = phi*y(1) - rho*y(1) + theta*y(3)*y(2) - alpha*y(3)*y(1);
    dydt(2) = eps*y(1) - theta*y(3)*y(2) - mu*y(2);
    dydt(3) = r*y(3)*(1 - b*y(3)) - gamma*y(1)*y(3);
end


