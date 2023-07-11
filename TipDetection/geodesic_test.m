clc;clear;close all;
x_size = 512; y_size = 511; z_size = 61;
% the last column is deleted becasue of the data property
% in the original data, the last column is almost empty
dbstop if error 
addpath C:\Users\mengf\Desktop\geodesic\geodesic_v3
% for_all bound_all node_num bound_num win_size

win_size = 15;
neighbor = cat(3, [0 0 0; 0 1 0; 0 0 0],...
                  [0 1 0; 1 1 1; 0 1 0],...
                  [0 0 0; 0 1 0; 0 0 0]);
% calculate distance
step = (win_size-1)/2;
% center = [(win_size+1)/2 (win_size+1)/2 (win_size+1)/2];
% center_idx = sub2ind([win_size win_size win_size],center(1),center(2),center(3));
              
load fore_all.mat
tic;
fore_all = double(fore_all);
bound_temp = fore_all - imerode(fore_all, neighbor);
bound_all = zeros(size(fore_all));
bound_all(step+1:end-step,step+1:end-step,step+1:end-step) = bound_temp(step+1:end-step,step+1:end-step,step+1:end-step);

fore_ind = find(fore_all);
node_num = length(fore_ind);
fore_all(fore_ind) = 1:node_num;
bound_ind = find(bound_all);
bound_num = length(bound_ind);
toc
% [x_ind,y_ind,z_ind] = ind2sub(size(bound_all),find(bound_all));
% bound_list = [x_ind y_ind z_ind];


distance_map = main(node_num, fore_all, bound_num, bound_all, [6 7], [1 1 1], 1);
toc

map_3 = max(distance_map,[],3);
% imagesc(map_3);colorbar;
imagesc(map_3,'AlphaData',~isnan(map_3));
set(gca,'color',0*[1 1 1]);
colorbar;
toc