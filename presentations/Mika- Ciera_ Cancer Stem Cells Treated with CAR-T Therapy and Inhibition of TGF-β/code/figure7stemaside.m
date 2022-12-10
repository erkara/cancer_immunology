%figure3a results
%HDLM-2 cell-line with CAR-T 123 therapy 
clc;clearvars;
dt = 0.01;
%firs injection, second challenge and max time
T_inject1 = 7;
T_max = 80;
tspan = 0:dt:T_max;
index1 = find(tspan == T_inject1);
%tumor grows without intervention
t1 = tspan(1:index1);
T0=3.4e7;
S0=0.01*T0;
E0=5e6;
C0=0;
B0=5;
R0=2.53;

[~,y] = ode23s(@(t,y) GetODE(t,y),t1,[T0,S0,E0,C0,B0,R0]);
%store all variables
T1 = y(:,1);
S1 = y(:,2);
E1 = y(:,3);
C1 = y(:,4);
B1 = y(:,5);
R1 = y(:,6);
% 
%t = 7 days, inject 2e6 CART-Cells  
t2 = tspan(index1+1:end);
T0=T1(end);
S0=S1(end);
E0=E1(end)+5e6;
C0=C1(end);
B0=B1(end);
R0=R1(end);

[~,y] = ode23s(@(t,y) GetODE(t,y),t2,[T0,S0,E0,C0,B0,R0]);
T = [T1;y(:,1)];
S = [S1;y(:,2)];
E = [E1;y(:,3)];
C = [C1;y(:,4)];
B = [B1;y(:,5)];
R = [R1;y(:,6)];


figure;
plot(tspan,T)
ylabel("Tumor Cells")
ylim([0,6e8])
saveas(gcf,'fgi7a1.jpg','jpg')


figure;
plot(tspan,S)
ylabel("Cancer Stem Cells")
ylim([0,2.8e6])
saveas(gcf,'fgi7a2.jpg','jpg')

figure;
plot(tspan,C+E)
xlabel("t")
ylabel("Tumor Effcetive cells")
ylim([0,6e6])
saveas(gcf,'fgi7a3.jpg','jpg')


figure;
plot(tspan,R)
xlabel("t")
ylabel("Regulatory T-cells")
ylim([0,2.5e5])
saveas(gcf,'fgi7a4.jpg','jpg')


figure;
plot(tspan,B)
xlabel("t")
ylabel("TGF-beta")
ylim([0,60])
saveas(gcf,'fgi7a5.jpg','jpg')


function dydt = GetODE(t,y)
    dydt = zeros(6,1);
    %parameters from the paper
    p = 0.1351;
    v0 = 0.3628;
    k = 1e9;
    muCS = 1e-6;
    eta = 0.0454;
    muET = 2e-5;
    c1 = 100;
    b = 0.4843;
    f = 0.0011;
    c3 = 0.1224;
    muRE = 1e-5;
    delE = 1e-5;
    r = 0.0125;
    muRC = 1e-5;
    delCT = 1e-4;
    a = 0.5917;
    c2 = 1.875e9;
    alphaRB = 1e-8;
    delB = 7e-4;
    delR = 1e-5;
    T = y(1);
    S = y(2);
    E = y(3);
    C = y(4);
    B = y(5);
    R = y(6);
    
   
    
    dydt(1) = eta*T*(1-((T+S)/k)) + (1-p)*v0*(1-((T+S)/k))*S - muET*((T*E)/(1+c1*B));
    dydt(2) = p*v0*(1 -((T+S)/k))*S - muCS*C*S;
    dydt(3) = b*E + f*E*T/(1+c3*T*B) - muRE*R*E - delE*E-r*E;
    dydt(4) = b*C + f*(C*S/(1+c3*T*B)) - muRC*R*C - delCT*C;
    dydt(5) = (a*T^2)/(c2+T^2) + alphaRB*R - delB*B;
    dydt(6) = r*E - delR*R;
end