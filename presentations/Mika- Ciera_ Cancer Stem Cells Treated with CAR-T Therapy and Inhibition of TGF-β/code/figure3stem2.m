clc;close all;clearvars;

% %firs injection, second challenge and max time
% T_inject1 = 42;
% T_inject2 = 250;
dt = 0.1; T_max = 500000; tspan = 0:dt:T_max;
% index1 = find(tspan == T_inject1);
% index2 = find(tspan == T_inject2);

T0=3.4e7;
S0=0.01*T0;
E0=0;
C0=0;
B0=5;
R0=2.5;

[~,y] = ode23s(@(t,y) GetODE(t,y),tspan,[T0,S0,E0,C0,B0,R0]);
%store all variables
T = y(:,1);
S = y(:,2);
E = y(:,3);
C = y(:,4);
B = y(:,5);
R = y(:,6);


% %run for 500 days
% figure;
% plot(tspan,T)
% ylabel("Tumor Cells")
% %ylim([0,1e9])
% xlim([0,T_max])
% saveas(gcf,'fgi3a.jpg','jpg')
% 
% %run for 500 days
% figure;
% plot(tspan,S)
% ylabel("Cancer Stem Cells")
% %ylim([0,12e6])
% xlim([0,T_max])
% saveas(gcf,'fgi3b.jpg','jpg')

% %run for 5,000 days
% figure;
% plot(tspan,B)
% ylabel("TGF-")
% ylim([0,1000])
% saveas(gcf,'fgi3c.jpg','jpg')

%run for 500,000 days
figure;
plot(tspan,R)
ylabel("Regulatory T-cells")
ylim([0,3])
saveas(gcf,'fgi3d.jpg','jpg')


% function dydt = GetHDLM(t,y)
%     dydt = zeros(6,1);
%     %parameters from the paper
%     p = 0.1351;
%     v0 = 0.3628;
%     k = 10e9;
%     muCS = 10e-6;
%     eta = 0.0454;
%     muET = 2*10e-5;
%     c1 = 100;
%     b = 0.4843;
%     f = 0.0011;
%     c3 = 0.1224;
%     muRE = 10e-5;
%     delE = 10e-5;
%     r = 0.0125;
%     muRC = 10e-5;
%     delCT = 10e-4;
%     a = 0.5917;
%     c2 = 1.875*10e9;
%     alphaRB = 10e-8;
%     delB = 7*10e-4;
%     delR = 10e-5;
%     
%     dydt(1) = p*v0*(1-(y(1)+y(2)/k))*y(2) - muCS*y(4)*y(2);
%     dydt(2) = eta*y(1)*(1-((y(1)+y(2))/k)) + (1-p)*v0*(1-(y(1)+y(2)/k))*y(2) - muET*(y(1)*y(3)/1+c1*y(5));
%     dydt(3) =  b*y(3)+(f*y(3)*y(1)/1+c3*y(1)*y(5))- muRE*y(6)*y(3)-delE*y(3)-r*y(3);
%     dydt(4) = b*y(4)+(f*(y(4)*y(2)/1+c3*y(1)*y(5)))- muRC*y(6)*y(4)-delCT*y(4);
%     dydt(5) = (a*y(1)^2/c2+y(1)^2)+alphaRB*y(6) - delB*y(5);
%     dydt(6) = r*y(3)-delR*y(6);
% end

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