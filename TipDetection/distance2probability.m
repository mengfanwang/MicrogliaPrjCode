clc;clear;close all;
x_size = 512; y_size = 511; z_size = 61;
% the last column is deleted becasue of the data property
% in the original data, the last column is almost empty
dbstop if error 

neighbor = cat(3, [0 0 0; 0 1 0; 0 0 0],...
                  [0 1 0; 1 1 1; 0 1 0],...
                  [0 0 0; 0 1 0; 0 0 0]);
              
step = 21;
              
for tt = 1:1
ind = num2str(1000+tt);
ind = ind(2:4);
fore_all = load(['D:\dropbox\foreground\' ind '.mat']);
fore_all = logical(fore_all.labelMap);
fore_all = fore_all(:,1:y_size,:);




% load fore_all.mat
% remove small regions
com = bwconncomp(fore_all);
temp = zeros(size(fore_all));
com_size = cellfun(@length,com.PixelIdxList);
for ii = 1:com.NumObjects
    if com_size(ii) > 2000 % small region threshold
        temp(com.PixelIdxList{ii}) = fore_all(com.PixelIdxList{ii});
    end
end
fore_all = temp;

fore_all = double(fore_all);
bound_temp = fore_all - imerode(fore_all, neighbor);
% bound_all = zeros(size(fore_all));
% bound_all(step+1:end-step,step+1:end-step,step+1:end-step) = bound_temp(step+1:end-step,step+1:end-step,step+1:end-step);

fore_all = padarray(fore_all,[step step step],'symmetric','both');
bound_all = zeros(size(fore_all));
bound_all(step+1:end-step,step+1:end-step,step+1:end-step) = bound_temp;
fore_ind = find(fore_all);
node_num = length(fore_ind);
fore_all(fore_ind) = 1:node_num;
bound_ind = find(bound_all);
bound_num = length(bound_ind);

tic;
distance_map = geodesic_v3(node_num, fore_all, bound_num, bound_all,1:step, [1 1 1], 2);
toc
distance_map = distance_map(step+1:end-step,step+1:end-step,step+1:end-step,:);
for ww = 1:step
    ground_truth = load(['D:\dropbox\distance_map\' num2str(ww) '\1.mat']);
    ground_truth = ground_truth.distance_map;
    sum(abs(distance_map(:,:,:,ww) - ground_truth),'all','omitnan')
end
end


% % p_value = 0.05;
% % distance = distance_map(distance_map<100&distance_map > 0.01);
% % distance = sort(distance,'descend');
% threshold = distance(floor(p_value*length(distance)));
% 
% result = distance_map >= threshold;
% im = zeros(x_size,y_size,3,z_size);
% for zz = 1:z_size
%     im(:,:,1,zz) = ~isnan(distance_map(:,:,zz)) * 128 + result(:,:,zz)*127;
%     im(:,:,2,zz) = ~isnan(distance_map(:,:,zz)) * 128;
%     im(:,:,3,zz) = ~isnan(distance_map(:,:,zz)) * 128;
% end
% im = uint8(im);
% map_3 = max(distance_map,[],3);
% imagesc(map_3,'AlphaData',~isnan(map_3));
% set(gca,'color',0*[1 1 1]);
% colorbar;

% load map_31.mat
% load fore_all.mat
% com = bwconncomp(fore_all);
% temp = nan(size(fore_all));
% com_size = cellfun(@length,com.PixelIdxList);
% for ii = 1:com.NumObjects
%     if com_size(ii) > 2000 % small region threshold
%         temp(com.PixelIdxList{ii}) = distance_map(com.PixelIdxList{ii});
%     end
% end
% % writeColormap(temp,'./temp_data/temp')
% distance_map = temp;
% x = distance_map(distance_map< 100 & distance_map > 0.01);

% Expontential fitting
% B = 40;
% x_bar = mean(b(b<B));
% lambda = truncatedExponential(x_bar,B);
% y = lambda*exp(-lambda*[0:64]);
% plot(y/sum(y)); hold on;

% % Gamma fitting
% B = min(8,max(x)*0.8);
% B = round(B*10)/10;
% x_bar = mean(x(x<B));
% logx_bar = mean(log(x(x<B)));
% fun = @(x) gammaDistributionTruncated(x,B,x_bar,logx_bar);
% x0 = [1 x_bar];
% result = fsolve(fun, x0)
% % plot histogram
% figure(1);
% h = histogram(x(x<B),'BinWidth',0.1);
% p = h.Values;
% h = histogram(x,'BinWidth',0.1);
% x_axis = 0:0.1:B;
% figure(2);
% bar(h.BinEdges(1:end-1)+0.05,h.Values/sum(p));hold on;
% cdf = gammainc(h.BinEdges/result(2),result(1))/gammainc(B/result(2),result(1));
% plot(h.BinEdges(1:end-1)+0.05, diff(cdf),'LineWidth',2);

function F = gammaDistributionTruncated(x, B, x_bar, logx_bar)
    a = x(1);
    b = x(2);
    gma = @(t,a) t.^(a-1).*exp(-t);
    gma_gradient = @(t,a) t.^(a-1).*exp(-t).*log(t);
    F(1) =  - a/b + x_bar/b^2 + (B/b)^a*exp(-B/b)/b/integral(@(t) gma(t,a),0,B/b);
    F(2) = logx_bar - log(b) - integral(@(t) gma_gradient(t,a),0,B/b)/integral(@(t) gma(t,a),0,B/b);
end

function lambda = truncatedExponential(x_bar,B)
lambda_left = 0;        
lambda_right = 1/x_bar;
lambda = (lambda_left + lambda_right)/2;
err = inf;
while err > 1e-8
    f = 1/lambda - B/(exp(lambda*B)-1) - x_bar;
    if f < 0
        lambda_right = lambda;
        lambda = (lambda_left + lambda_right)/2;
    elseif f > 0
        lambda_left = lambda;
        lambda = (lambda_left + lambda_right)/2;
    else
        break;
    end
    err = lambda_right - lambda_left;
end
end

function writeColormap(distance_map, filename)
distance_map(distance_map>8) = 8;
distance_map = ceil(distance_map*255/max(distance_map(:))-1e-8) + 1;
color = colormap;
[x,y,z] = size(distance_map);
map = zeros(x,y,3,z);
for ii = 1:x
    for jj = 1:y
        for kk = 1:z
            if ~isnan(distance_map(ii,jj,kk))
                map(ii,jj,1,kk) = color(distance_map(ii,jj,kk),1);
                map(ii,jj,2,kk) = color(distance_map(ii,jj,kk),2);
                map(ii,jj,3,kk) = color(distance_map(ii,jj,kk),3);
            end
        end
    end
end
tifwrite(uint8(map*255),filename);
end

