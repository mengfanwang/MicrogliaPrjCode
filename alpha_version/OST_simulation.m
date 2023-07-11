clc;clear;close all;

% mu2000 = dlmread('C:\Users\Mengfan Wang\Dropbox\code\Modify Series Data\mu2000.txt',',');
% sigma2000 = dlmread('C:\Users\Mengfan Wang\Dropbox\code\Modify Series Data\sigma2000.txt',',');

%%
% iter = 1000;
% m = 1000; n = 1000;
% value_label = zeros(m+n,2);
% data = normrnd(0,1,m+n,1);
% data = sort(data);
% data_small = data(1:2:end-1);
% data_large = data(2:2:end);
% value_label(1:n,1) = data_small;
% value_label(n+1:n+m,1) = data_large;
% value_label(n+1:n+m,2) = 1;
% value_label = sortrows(value_label, 1);
% label_large = logical(value_label(:,2));
% label_small = ~label_large;
%%
iter = 1000;
% m = 1000; n = 800;
l = zeros(iter,1);
label_large = rand(2000,1) < 0.5;
label_small = ~label_large;
m = sum(label_large);
n = sum(label_small);
for ii = 1:iter
    data = normrnd(0,1,m+n,1);
    data = sort(data);
    data_small = data(label_small);
    data_large = data(label_large);
    l(ii) = mean(data_large) - mean(data_small);
end
% fprintf('Text        : mu %f, sigma %f\n', mu2000(m,n), sigma2000(m,n));
fprintf('Simulation  : mu %f, sigma %f\n', mean(l), sqrt(var(l)));
[mu, sigma] = conventionalIntegralv1(data_large, data_small);
fprintf('Intergral   : mu %f, sigma %f\n', mu, sigma);
[mu, sigma] = ksegments_orderstatistics_v2(data_large, data_small);
fprintf('Acceleration: mu %f, sigma %f\n', mu, sigma);
[mu, sigma] = ordStatApproxKsec(data_large, data_small);
fprintf('Congchao    : mu %f, sigma %f\n', mu, sigma);

% pdfz = normpdf(norminv(n/(m+n),0,1));
% cdfz = normcdf(norminv(n/(m+n),0,1));
% mu = pdfz/(1-cdfz) + pdfz/cdfz;
% sigma = sqrt((1+norminv(n/(m+n),0,1)*pdfz/(1-cdfz)-(pdfz/(1-cdfz))^2)/m + (1-norminv(n/(m+n),0,1)*pdfz/cdfz-(pdfz/cdfz)^2)/n);
% fprintf('Zhengwei    : mu %f, sigma %f\n\n', mu, sigma);