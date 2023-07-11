 clc;clear;close all;
dbstop if error

%% truncated exponential distribution estimater
% lambda_initial = rand*100;
% x = exprnd(1/lambda_initial,100000,1);
% B = 12;
% x_bar = mean(x(x<B));
% 
% lambda_left = 0;
% lambda_right = 1/x_bar;
% lambda = (lambda_left + lambda_right)/2;
% err = inf;
% while err > 1e-8
%     f = 1/lambda - B/(exp(lambda*B)-1) - x_bar;
%     if f < 0
%         lambda_right = lambda;
%         lambda = (lambda_left + lambda_right)/2;
%     elseif f > 0
%         lambda_left = lambda;
%         lambda = (lambda_left + lambda_right)/2;
%     else
%         break;
%     end
%     err = lambda_right - lambda_left;
% end
% lambda_initial
% lambda

% x = 0:0.1:10;
% plot(x,log(x)-psi(x))
% gma = @(t,a) t.^(a-1).*exp(-t);
% gma_gradient = @(t,a) t.^(a-1).*exp(-t).*log(t);
% y = zeros(length(x));
% B = 6;
% for B = 1:5:11
% for ii = 1:length(x)
%     y(ii) = log(x(ii)) - integral(@(t) gma_gradient(t,x(ii)),0,B,'RelTol',1e-8,'AbsTol',1e-13)/...
%                          integral(@(t) gma(t,x(ii)),0,B,'RelTol',1e-8,'AbsTol',1e-13);
% end
% plot(x,y);hold on
% end

%% truncated gamma distribution estimater
% a_initial = 2;
% b_initial = 2;
% x = gamrnd(a_initial,b_initial,1000000,1);
load distance_map_max.mat
x = x*10;
B = 7;
x_bar = mean(x(x<B));
logx_bar = mean(log(x(x<B)));
fun = @(x) gammaDistributionTruncated(x,B,x_bar,logx_bar);
x0 = [2 x_bar/2];
result = fsolve(fun, x0)
% plot histogram
h = histogram(x(x<B),'BinWidth',0.1); hold on;
x = gamrnd(result(2),result(1),140000,1);
histogram(x,'BinWidth',0.1); hold on;
x_axis = 0:0.1:B;
figure(2);plot(x_axis(1:end-1)+0.05,h.Values/sum(h.Values));hold on;
cdf = gammainc(x_axis/result(2),result(1))/gammainc(B/result(2),result(1));
plot(x_axis(1:end-1)+0.05, diff(cdf));

%% truncated beta distribution estimator
clc;clear;close all
% a_initial = 10; b_initial = 30;
% x = betarnd(a_initial,b_initial,1000000,1);
load distance_map_max.mat
B = 0.66;
histogram(x(x<B));hold on;
logx_bar = mean(log(x(x<B)));
log1_x_bar = mean(log(1-x(x<B)));
% betaDistributionTrunacated([a_initial b_initial],B,logx_bar,log1_x_bar)
fun = @(x) betaDistributionTrunacated(x,B,logx_bar,log1_x_bar);
result = fsolve(fun, [3 5])
x = betarnd(result(1),result(2),139926,1);
histogram(x)
B = 0.34;
logx_bar = mean(log(x(x<B)))
log1_x_bar = mean(log(1-x(x<B)));

%% laplace distribution estimator
clc;clear;close all;
laplacernd = @(x,mu,b)  mu - b.*((x>0.5)*2-1).*log(1-2*abs(x-0.5));
% mu_initial = 10; b_initial = 10;
% x = rand(1000000,1);
% x = laplacernd(x,mu_initial,b_initial);
load distance_map_max.mat
% x = x(x>0.872);
% histogram(x,0.87:0.002:1);hold on;
% mu = 0.872;
% b = mean(abs(x-mu))
% y = rand(20000,1);
% y = laplacernd(y,mu,b);
% histogram(y,0.87:0.002:1)
h = histogram(x,0:0.01:1);hold on;
binValue = h.Values;
binValue(76:87) = binValue(76:87) - binValue(100:-1:89);
binValue(88:100) = 0;
bar(0.005:0.01:0.995, binValue);


function d = sampleBall
u = rand(1,3);
v = rand(1,3);
theta = 2*pi.*u;
phi = acos(2*v - 1);
x = sin(theta).*sin(phi);
y = cos(theta).*sin(phi);
z = cos(phi);
d = openGJK_v2([x; y; z]);
end




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

function F = gammaDistributionTruncated(x, B, x_bar, logx_bar)
    a = x(1);
    b = x(2);
    gma = @(t,a) t.^(a-1).*exp(-t);
    gma_gradient = @(t,a) t.^(a-1).*exp(-t).*log(t);
    F(1) =  - a/b + x_bar/b^2 + (B/b)^a*exp(-B/b)/b/integral(@(t) gma(t,a),0,B/b);
    F(2) = logx_bar - log(b) - integral(@(t) gma_gradient(t,a),0,B/b)/integral(@(t) gma(t,a),0,B/b);
end

function F = gammaDistribution(x,x_bar,logx_bar)
    a = x(1);
    b = x(2);
    gma = @(t,a) t.^(a-1).*exp(-t);
    gma_gradient = @(t,a) t.^(a-1).*exp(-t).*log(t);
    F(1) =  - a/b + x_bar/b^2;
    F(2) = logx_bar - log(b) - integral(@(t) gma_gradient(t,a),0,inf)/integral(@(t) gma(t,a),0,inf);
end
