clc;clear;close all;
% Dr. Wang and CT's method

load distance_map_max.mat

% trauncated fitting
% B = 0.66;
% histogram(x(x<B));hold on;
% logx_bar = mean(log(x(x<B)));
% log1_x_bar = mean(log(1-x(x<B)));
% fun = @(x) betaDistributionTrunacated(x,B,logx_bar,log1_x_bar);
% result = fsolve(fun, [3 5])
% x = betarnd(result(1),result(2),139926,1);
% histogram(x)

logx_bar = mean(log(x));
log1_x_bar = mean(log(1-x));
fun = @(x) betaDistribution(x,logx_bar,log1_x_bar);
result = fsolve(fun, [3 5]);
a = result(1); b = result(2);
% histogram(x);hold on;
% x = betarnd(a,b,length(x),1);
% histogram(x);

threshold = betainv(0.95,a,b);
p_value = sum(x>threshold)/length(x)

for ii = 1:99
    alpha = ii/100;
    threshold = betainv(alpha,a,b);
    p_value(ii) = sum(x<threshold)/length(x);
end
figure(2);
plot(p_value);hold on;
plot(0.01:0.01:0.99);

function F = betaDistributionTrunacated(x,B,logx_bar,log1_x_bar)
    a = x(1);
    b = x(2);
    
    beta_g2a = @(t,a,b) t.^(a-1).*(1-t).^(b-1).*log(t);
    beta_g2b = @(t,a,b) t.^(a-1).*(1-t).^(b-1).*log(1-t);
    F(1) = logx_bar - integral(@(t) beta_g2a(t,a,b),0,B)/betainc(B,a,b)/beta(a,b);
    F(2) = log1_x_bar - integral(@(t) beta_g2b(t,a,b),0,B)/betainc(B,a,b)/beta(a,b);
end


function F = betaDistribution(x,logx_bar,log1_x_bar)
    a = x(1);
    b = x(2);
    F(1) = logx_bar - psi(a) + psi(a+b);
    F(2) = log1_x_bar - psi(b) + psi(a+b);
end