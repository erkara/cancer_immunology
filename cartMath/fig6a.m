%plot fig6a, percentage of the virtual cohort surviving less than the given
%number of days. First run vp_generate code 
clc;clearvars;format bank;
data = table2array(readtable("results/vp.csv"));
sur_days = sort(data(:,10));
d = 1:300;
len = length(sur_days);
percents = zeros(1,length(d));
for i=1:length(d)
    percents(i) = round(100*(len-length(find(sur_days > d(i))))/len,2); 
end
plot(d,flip(percents),'r-',linewidth=1)
xlabel('percent survival');
ylabel('days')
grid on;
%check for different doses


