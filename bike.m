% This is the Matlab code used for the real data analysis in
% "AN RKHS APPROACH FOR PIVOTAL INFERENCE IN FUNCTIONAL LINEAR REGRESSION"
% by Holger Dette and Jiajun Tang
%
% The data analyzed is the Bike-sharing data of Captial Bike Sharing (CBS), obtained from R package ISLR2.
% References for the Bike-sharing data:
% Fanaee-T, H. and Gama, J. (2014). Event labeling combining ensemble detectors and background knowledge. Progress in Artificial Intelligence, 2, 113–127.
% James, G., Witten, D., Hastie, T. and Tibshirani, R. (2021). An Introduction to Statistical Learning: With Applications in R. New York: springer.
%
% Please download and source the Matlab package Chebfun available at https://www.chebfun.org/



clear
% Get data
% Bike
bikeshare = readtable('bikershare.txt');
bikeshare.Properties.VariableNames={'season', 'mnth', 'day',...
    'hr', 'holiday', 'weekday', 'workingday', 'weathersit', 'temp', 'atemp',...
    'hum', 'windspeed', 'casual', 'registered', 'bikers'};

bikeshare= bikeshare(:,{'day','hr','workingday','temp','windspeed','bikers'});
bikeshare = table2array(bikeshare);
workday_id = find(bikeshare(:,3)==1);
workday = unique(bikeshare(workday_id,1));
bikeshare = bikeshare(workday_id,:);
bikes = zeros(1,250);
t = zeros(100,250);

t = 0:0.01:1;
temp = zeros(101,250);
wind = zeros(101,250);

hr0 = cell(1,250);
temp0 = cell(1,250);
wind0 = cell(1,250);

for i = 1:250
    day = workday(i);
    id = (bikeshare(:,1)==day);
    bikes(i) = sum(bikeshare(id,6));
    hr0{i} = bikeshare(id,2)/23;    
    temp0{i} = bikeshare(id,4);
    wind0{i} = bikeshare(id,5);
    f = fit(hr0{i},temp0{i}, 'fourier3');
    temp(:,i)=f(t);
    f = fit(hr0{i},wind0{i}, 'fourier3');
    wind(:,i)=f(t);
end


% Compute the quantiles of the pivotal statistic W
Q = 25;
nu0 = 1/2;
dt = (1-nu0) / Q;
nn = 10000;
ww = zeros(nn, 1);
nuq = nu0 + (1-nu0)/Q*(1:Q)';
for i = 1:nn
    b = zeros(Q, 1);
    b(1) = sqrt(nu0+dt) * randn;
    for j = 2:Q
        b(j) = b(j-1)+sqrt(dt)*randn;
    end
    %     temp = nuq .* (nuq .* b - b(end));
    t = nuq .* b - b(end);
    t = norm(t)^2;
    t = sqrt(t/Q*(1-nu0));
    ww(i) = b(end) / t;
end
Wq = quantile(ww, [0.99, 0.95, 0.9]);
% clear ww

bikes = (bikes-min(bikes))/(max(bikes)-min(bikes));




%% Temperature
Delta_temp = 1:0.01:1.7;
rj_temp = zeros(3,length(Delta_temp));
gridx = 0:0.01:1;
[T,V] = TandV(temp',bikes',gridx,Q,nu0);
% Temperature
for a = 1:3
    for d = 1:length(Delta_temp)
        rj_temp(a,d) = (T>Wq(a)*V+Delta_temp(d));
    end
end

% Confidence intervals alpha=0.90
% One sided 
[0, T+quantile(ww, [0.90])*V]
% Two sided
[max(0,T-quantile(ww, [0.95])*V), T+quantile(ww, [0.95])*V]

% Confidence intervals alpha=0.95
% One sided 
[0, T+quantile(ww, [0.95])*V]
% Two sided
[max(0,T-quantile(ww, [0.975])*V), T+quantile(ww, [0.975])*V]

% Confidence intervals alpha=0.95
% One sided 
[0, T+quantile(ww, [0.99])*V]
% Two sided
[max(0,T-quantile(ww, [0.995])*V), T+quantile(ww, [0.995])*V]


%% Wind speed
Delta_wind = 9.7:0.1:11.4;
rj_wind = zeros(3,length(Delta_wind));
 [T,V] = TandV(wind',bikes',gridx,Q,nu0);
for a = 1:3
    for d = 1:length(Delta_wind)
        rj_wind(a,d) = (T>Wq(a)*V+Delta_wind(d));
    end
end

% Confidence intervals alpha=0.90
% One sided 
[0, T+quantile(ww, [0.90])*V]
% Two sided
[max(0,T-quantile(ww, [0.95])*V), T+quantile(ww, [0.95])*V]

% Confidence intervals alpha=0.95
% One sided 
[0, T+quantile(ww, [0.95])*V]
% Two sided
[max(0,T-quantile(ww, [0.975])*V), T+quantile(ww, [0.975])*V]

% Confidence intervals alpha=0.95
% One sided 
[0, T+quantile(ww, [0.99])*V]
% Two sided
[max(0,T-quantile(ww, [0.995])*V), T+quantile(ww, [0.995])*V]
