%things to investigate
%1- try different challange amounts
%2- try different challange days
%3- initial tumor burden T0

clc;close all;clearvars;
%solve your odes here
dt=0.1;Tmax = 40;tspan = 0:dt:Tmax;
%injecttion params
%challange with new tumor cells
dayinject = 10;
index = find(tspan==dayinject);
challenge = 2e6;
type = "nn";
%solve ode until day 10
t1 = tspan(1:index);
T0 = 1e4/2;    
N0 = 1e2;
L0 = 1e1;

[~,y]=ode23s(@(t,y) GetODEHuman(t,y,type),t1,[T0,N0,L0]);
T = y(:,1);
N = y(:,2);
L = y(:,3);

%inejct challenge number of cells and solve ode after day10
t2 = tspan(index+1:end);
T0 = T(end) + challenge;    
N0 = N(end);
L0 = L(end);
[~,y]=ode23s(@(t,y) GetODEHuman(t,y,type),t2,[T0,N0,L0]);
T = [T ; y(:,1)];
N = [N ; y(:,2)];
L = [L ; y(:,3)];


semilogy(tspan,T,'r-',LineWidth=2)
hold on
semilogy(tspan,N,'b-.',LineWidth=2)
hold on
semilogy(tspan,L,'g--',LineWidth=2)
ylim([1e2,1e9]);
legend("Tumor","NKcells","CD8+cells")
grid on;
xlabel ('Rechallenge Day 10')
ylabel ( 'Total Cell Populatiton')


function dydt = GetODEHuman(t,y,type)
    dydt = zeros(3,1);%[T,N,L]
    a = 5.14e-1; 
    b = 1.02e-9;
    if type=="nl"
        c = 3.23e-7;
        d = 1.43;
        lambda = 5.8e-1;
        s = 2.73;
        g = 2.5e-2;
        j = 3.75e-2;
    elseif type=="nn"
        c = 3.23e-7;
        d = 3.6;
        lambda = 4.6e-1;
        s = 1.61;
        g = 2e-1;
        j = 3.75e-2;
    else
        disp("enter a valid ligand type")
    end
    
    sigma = 1.3e4;
    f = 4.12e-2; 
    h = 2.02e7; 
    k = 2.02e7;
    m = 2.00e-2; 
    q = 3.42e-10; 
    p = 1.0e-7; 
    r = 1.1e-7;

    %ODE paramaters
    T = y(1);
    N = y(2);
    L = y(3);
    
    %days
    r1 = (L/T)^lambda;
    D = d*(r1/(s+r1))*T;
    
    dydt(1) = a*T*(1-b*T) - c*N*T - D; %T
    
    dydt(2) = sigma - f*N + ((g*T^2)/(h+T^2))*N - p*N*T; %N
    
    dydt(3) = -m*L + ((j*D^2)/(k+D^2))*L - q*L*T + r*N*T;%L

end

