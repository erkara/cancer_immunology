%suggested solution of Q3 in hw1
clear all;close all;clc;
%read the data
data = table2array(readtable('tumor.csv'));
days = data(:,1);
tumorsize = data(:,2);
%model names
models = ["Exponential","Mendelsohn","Logistic",...
    "Linear","Surface","Gompertz","Bertalanffy"];
params = [1,2,2,2,2,3,2];
%solve from the first to the last day of experiments
tspan = days(1):1:days(end);
%this is the first tumor size in your dataset,change it accordingly
y0 = 225;
%let's store the solutions in the following matrix;
sol = zeros(length(tspan),length(models));

%loop over model names, get and solve ode and plot one at a time
for i=1:length(models)
    %solve the model
    [t,y] = ode45(@(t,y) GetODE(t,y,models(i)),tspan,y0);
    %register the solutions to sol
    sol(:,i) = y;
    %using the first 7 data points
    n = 7;
    %number of parameters in the related model
    K = params(i);
    %compute ssr
    ssr = sum((y(days(1:n))-tumorsize(1:n)).^2);
    %compute AIC
    aic = n*log(ssr/n)+(2*(K+1)*n)/(n-K-2);
    fprintf("model: %s ssr: %0.0f aic: %0.0f\n",models(i),ssr,aic);
    %plot the 
    plot(t,y)
    ylim([0,1600]);
    xlim([0,65]);
    hold on;
end

grid on;
plot(days(1:7),tumorsize(1:7),'ks','Markersize',10,'MarkerFaceColor','k')
title('Figure-1');
xlabel('time');
ylabel('tumor size');
legend([models,'data'],Location="northwest");


%plot the solutions for all days
figure();
for i=1:length(models)
    plot(tspan,sol(:,i))
    hold on
end
plot(days,tumorsize,'ks','Markersize',10,'MarkerFaceColor','k');
legend([models,'data'],Location="northwest");
grid on;
title('Figure-2');
xlabel('time');
ylabel('tumor size');




function dydt = GetODE(t, y, model)
        if model=="Exponential"
            a= 0.0262;
            dydt = a*y;
        elseif model=="Mendelsohn"
            a = 0.286;
            b = 0.616;
            dydt = a*y^b;
        elseif model=="Logistic"
            a = 0.0370;
            b = 2000;
            dydt = a*y*(1-(y/b));
        elseif model=="Linear"
            a = 58.7;
            b = 1690;
            dydt = (a*y)/(y+b);
        elseif model=="Surface"
            a = 0.265;
            b = 506;
            dydt = (a*y)/(y+b)^(1/3);
        elseif model=="Gompertz"
            a = 0.279;
            b = 13900;
            c = 12000;
            dydt = a*y*log(b/(y+c));
        elseif model=="Bertalanffy"
            a = 0.306;
            b = 0.0119;
            dydt = a*y^(2/3)-b*y;
        end
end
