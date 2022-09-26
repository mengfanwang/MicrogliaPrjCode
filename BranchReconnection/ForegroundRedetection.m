clc;clear;close all;

% it will re-detect foreground based on the order-statistic thresholding
% result, connect broken branches and remove connected tips

% Algorithm:
% Input: foreground (from order-stattistic thresholding or sth else)
% 1. Choose a threshold from top-to-down thresholds and find candidate foreground
% 2. If a candidate component connect multiple components:
%       find the minimum spanning tree of them
%       For every connected nodes:
%           find the path by region growing
%           use t-test to accept the path or not

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';
data_all = load([file_name 'stabilization_c2\stabilization_data.mat']);
data_all = data_all.data;
[x,y,z,t] = size(data_all);

conn = 8;
gaussian_sigma = 0.8;

tt = 1;
ind = num2str(1000+tt);
ind = ind(2:4);
fore_all = load([file_name 'foreground_c2\' ind '.mat']);
fore_all = fore_all.foreground;

p_map = zeros(x,y,z);
tic;
for zz = 34:34

data = data_all(:,:,zz,tt);
data_smoothed = imgaussfilt(data,gaussian_sigma);
foreground = fore_all(:,:,zz);
% imshow(foreground);
com_before = bwconncomp(foreground,conn);
p_map_temp = zeros(x,y);


for thre = 20:-1:1 
    binaryData = zeros(x,y);
    binaryData(data_smoothed>thre) = 1;    
    binaryData = binaryData .* (1 - foreground);
%     imshow(binaryData);
    com_candidate = bwconncomp(binaryData,conn);
    for ii = 1:com_candidate.NumObjects
        foreground_temp = foreground;
        foreground_temp(com_candidate.PixelIdxList{ii}) = 1;
        com_after = bwconncomp(foreground_temp,conn);
        if com_after.NumObjects < com_before.NumObjects
            % get connected foregrounds
            com_list = [];
            for jj = 1:com_before.NumObjects
                foreground_temp = zeros(x,y);
                foreground_temp(com_before.PixelIdxList{jj}) = 1;
                foreground_temp(com_candidate.PixelIdxList{ii}) = 1;
                com_temp = bwconncomp(foreground_temp,conn);
                if com_temp.NumObjects == 1
                    com_list = [com_list jj];
                end
            end
%             % debug
% %             figure(2);imshow(foreground_temp);
%             com_list
%             fore_debug = zeros(x,y);
%             for dd = 1:length(com_list)
%                 fore_debug(com_before.PixelIdxList{com_list(dd)}) = com_list(dd);
%             end
%             ii
%             % debug off
            % generate min span tree
            mask = zeros(size(data));
            mask(com_candidate.PixelIdxList{ii}) = 1;
            ms_tree = genMinspantree(data, mask, com_before, com_list);
            for jj = 1:size(ms_tree.Edges.EndNodes,1)
                com1 = str2num(ms_tree.Edges.EndNodes{jj,1});
                com2 = str2num(ms_tree.Edges.EndNodes{jj,2});
                [pvalue, mask_final] = candidateTest(data, mask, com_before, com1, com2);
                fprintf('%d %d %f\n',com1,com2,pvalue);
                if pvalue > 0.05
%                     p_map(mask_final,zz) = pvalue;
                    p_map_temp(mask_final) = max(p_map_temp(mask_final), pvalue);
                end
            end
        end
    end
end
p_map(:,:,zz) = p_map_temp;
zz
toc
end
% tifwrite( (p_map>0)+fore_all, 'temp')

function [pvalue, mask3] = candidateTest(data,mask,com_before,com1,com2)
%     com1 = 16;com2 = 20;
    % find shortest path between two components
    cost_map = 1./(data+1);
    cost_map(~mask) = 1;
    mask1 = false(size(data));
    mask1(com_before.PixelIdxList{com1}) = 1;
    dist1 = graydist(cost_map, mask1);
    mask2 = false(size(data));
    mask2(com_before.PixelIdxList{com2}) = 1;
    dist2 = graydist(cost_map, mask2);
    dist_sum = dist1 + dist2;
    dist_sum = dist_sum - 1;
    dist_sum(logical(1-mask)) = inf;
    % region growing 
    mask3 = dist_sum < (min(dist_sum(:))*1.05);
%     mask3 = dist_sum == min(dist_sum(:));
    dist3 = bwdist(mask3, 'chessboard');
    mask4 = (dist3 < 3) & (mask1 | mask2);
%     [~, pvalue] = ttest2(data(mask3),data(mask4));
    pvalue = order_stat_test(sum(mask4(:)), sum(mask3(:)), mean(data(mask4)) - mean(data(mask3)));
%     imshow(mask1+mask2+mask3)
end

function t = genMinspantree(data, mask, com_before, com_list)
    cost_map = 1./(data+1);
    cost_map(~mask) = 1;
    g = graph();
    for ii = 1:length(com_list) - 1
        mask = zeros(size(data));
        mask(com_before.PixelIdxList{com_list(ii)}) = 1;
        dist_map = graydist(cost_map, logical(mask));
        for jj = ii+1:length(com_list)
            mask = zeros(size(data));
            mask(com_before.PixelIdxList{com_list(jj)}) = 1;
            g = addedge(g,num2str(com_list(ii)),num2str(com_list(jj)),min(dist_map(logical(mask))));
        end
    end
    t = minspantree(g);
end

function pvalue = order_stat_test(M,N,L)
    pdfz = normpdf(norminv(N/(M+N),0,1));
    cdfz = normcdf(norminv(N/(M+N),0,1));
    mu = pdfz/(1-cdfz) + pdfz/cdfz;
    sigma = sqrt((1+norminv(N/(M+N),0,1)*pdfz/(1-cdfz)-(pdfz/(1-cdfz))^2)/M + (1-norminv(N/(M+N),0,1)*pdfz/cdfz-(pdfz/cdfz)^2)/N);
    z_score = (L/sqrt(767) - mu)/sigma;  %variance for certian data
    pvalue = normcdf(z_score);
end

