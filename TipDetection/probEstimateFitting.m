clc;clear;close all;


tt = 1;
for step = 6:21
step
load(['D:/dropbox/distance_map/' num2str(step) '/' num2str(tt)]);

% x = distance_map(distance_map<100 & distance_map > 0.01);
% x(x>step-0.01) = step-0.01;
% x = x/step;

distance_map(distance_map < 0.01) = 0;
distance_map(distance_map>step-0.01) = step-0.01;
distance_map = distance_map/step;
if step == 6
    distance_map_max = distance_map;
    distance_map_arg = distance_map;
    distance_map_arg(~isnan(distance_map_arg)) = 6;
else
    distance_map_max = max(distance_map_max,distance_map);
    distance_map_arg(distance_map_max == distance_map) = step;
end
end
x = distance_map_max(~isnan(distance_map_max));
x = x(x>0);
x(x>0.89) = 0.89;
x = x/max(x);
figure(1);
h = histogram(x, 0:0.01:1);
p = h.Values/sum(h.Values);
% x = distance_map_arg(~isnan(distance_map_arg));
% figure(2);
% histogram(x(:));


% subplot(4,4,step-5);
% histogram(x, 0.9:0.01:1);
% xlabel(step);

% B = max(x)*0.7;
% B = ceil(B*10)/10;
% x_bar = mean(x(x<B));
% logx_bar = mean(log(x(x<B)));
% fun = @(x) gammaDistributionTruncated(x,B,x_bar,logx_bar);
% x0 = [1 x_bar];
% result = fsolve(fun, x0)
% % plot histogram
% figure(1);
% h = histogram(x(x<B),'BinWidth',0.1);
% p1 = h.Values;
% h = histogram(x,'BinWidth',0.1);
% x_axis = 0:0.1:B;
% figure(2);subplot(3,3,step-9);
% bar(h.BinEdges(1:end-1)+0.05,h.Values/sum(p1));hold on;xlabel(2*step+1);
% cdf = gammainc(x_axis/result(2),result(1))/gammainc(B/result(2),result(1));
% plot(x_axis(1:end-1)+0.05, diff(cdf),'LineWidth',2);
% 
% x = step - x;
% B = max(x)*0.3;
% B = ceil(B*10)/10;
% x_bar = mean(x(x<B));
% logx_bar = mean(log(x(x<B)));
% fun = @(x) gammaDistributionTruncated(x,B,x_bar,logx_bar);
% x0 = [7 x_bar/7];
% result = fsolve(fun, x0)
% % plot histogram
% figure(1);
% h = histogram(x(x<B),'BinWidth',0.1);
% p2 = h.Values;
% h = histogram(x,'BinWidth',0.1);
% x_axis = 0:0.1:B;
% figure(2);subplot(3,3,step-9);
% cdf = gammainc(x_axis/result(2),result(1))/gammainc(B/result(2),result(1));
% plot(step-x_axis(1:end-1)-0.05, diff(cdf)*sum(p2)/sum(p1),'LineWidth',2);


function F = gammaDistributionTruncated(x, B, x_bar, logx_bar)
    a = x(1);
    b = x(2);
    gma = @(t,a) t.^(a-1).*exp(-t);
    gma_gradient = @(t,a) t.^(a-1).*exp(-t).*log(t);
    F(1) =  - a/b + x_bar/b^2 + (B/b)^a*exp(-B/b)/b/integral(@(t) gma(t,a),0,B/b);
    F(2) = logx_bar - log(b) - integral(@(t) gma_gradient(t,a),0,B/b)/integral(@(t) gma(t,a),0,B/b);
end


% population = 1:125;
% population(63) = [];
% b = zeros(85000,1);
% for bb = 1:85000
% choice = randsample(population,4,true);
% choice = unique(choice);
% [x_ind,y_ind,z_ind] = ind2sub([5 5 5],choice);
% line_ind =[x_ind; y_ind; z_ind;];
% b(bb) = openGJK_v2(line_ind - [3; 3; 3;]);
% end
% b = b(b<100 & b > 0.01);
% histogram(b,'BinWidth',0.01);