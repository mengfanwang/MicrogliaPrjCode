clc;clear;
% get local maximum and component tree
addpath D:\MatlabTools\
addpath ..\
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
foreground_name = fullfile(file_name, 'foreground_iter1_jointThresholding3\0.05_5');

%% parameter setting
gaussian_sigma = 1;
max_dist = 15;
min_size = 20;
priCvt_thre = 0.5;
conn = ones(3,3,3);
%%
% get data          
load([file_name '\data_c2']);
load([file_name '\stabilizeFunction_c2']);
data_all = double(data2);
[x,y,z,t] = size(data_all);
std = sqrt(max(variance)/3);
data_all = interp1(0:255,stabilizeFunction,data_all);

data_all_new = zeros(size(data_all));
for tt = 1:t
    [x_mesh,y_mesh,z_mesh] = meshgrid(-4:x+5,-4:y+5,-4:z+5);
    [x_q,y_q,z_q] = meshgrid(1:x,1:y,1:z);
    if tt > 1
        load([file_name '\result_gradientDescent\result_gradientDescent_' num2str(tt-1) '_0.01_0.2_1500.mat']);
        data1 = data_all(:,:,:,tt-1);
        data1_pad = padarray(data1,[5 5 5],'replicate'); 
        data1_forward = interp3(x_mesh,y_mesh,z_mesh,data1_pad,x_q+uy,y_q+ux,z_q+uz);
    end
    if tt < t
        load([file_name '\result_gradientDescent_inverse\result_gradientDescent_' num2str(tt+1) '_0.01_0.2_1500.mat']);
        data1 = data_all(:,:,:,tt+1);
        data1_pad = padarray(data1,[5 5 5],'replicate'); 
        data1_backward = interp3(x_mesh,y_mesh,z_mesh,data1_pad,x_q+uy,y_q+ux,z_q+uz);
    end
    if tt == 1
        data_all_new(:,:,:,tt) = (data_all(:,:,:,tt) + data1_backward)/2;
    elseif tt == t
        data_all_new(:,:,:,tt) = (data_all(:,:,:,tt) + data1_forward)/2;
    else
        data_all_new(:,:,:,tt) = (data_all(:,:,:,tt) + data1_forward + data1_backward)/3;
    end
end
data_all = data_all_new;

%%  read data, foreground and build cost map
for tt = 2:t-1
    fprintf('Processing %d/%d...\n', tt, t);
data = data_all(:,:,:,tt);
% load foreground and remove small objects
foreground = tifread(fullfile(foreground_name, [num2str(tt) '.tif']));
foreground = logical(foreground);
com = bwconncomp(foreground,conn);
for ii = 1:com.NumObjects
    if length(com.PixelIdxList{ii}) < min_size
        foreground(com.PixelIdxList{ii}) = 0;
    end
end
% get principal curvature
eig_all = principalCv3d(data, gaussian_sigma);
eig_all(eig_all < 0) = 0;
cost_map = 1./eig_all.^4;
cost_map(cost_map > 1e8) = 1e8;
cost_map(foreground) = 1e8;

%% find the shortest path
g = graph;
com = bwconncomp(foreground,conn);
fore_candidate = cell(1,1);
for ii = 1:com.NumObjects
    mask_temp = zeros(x,y,z);
    mask_temp(com.PixelIdxList{ii}) = 1;
    dist_temp = bwdist(mask_temp);
%     fprintf('%d/%d ', ii, com.NumObjects);
    for jj = ii+1:com.NumObjects
        if min(dist_temp(com.PixelIdxList{jj})) < max_dist
            % get the range of cells
            [x_min, y_min, z_min, x_max, y_max, z_max] = getUnionRange(foreground, com, ii,jj);
            foreground_temp = foreground(x_min:x_max, y_min:y_max, z_min:z_max);
            % find shortest path
            % mask1&2: two foreground %mask3: shortest path
            mask1 = false(x,y,z);
            mask1(com.PixelIdxList{ii}) = 1;
            mask2 = false(x,y,z);
            mask2(com.PixelIdxList{jj}) = 1;
            mask1 = mask1(x_min:x_max, y_min:y_max, z_min:z_max);
            mask2 = mask2(x_min:x_max, y_min:y_max, z_min:z_max);
            dist1 = graydist(cost_map(x_min:x_max, y_min:y_max, z_min:z_max), mask1, 'quasi-euclidean');
            dist2 = graydist(cost_map(x_min:x_max, y_min:y_max, z_min:z_max), mask2, 'quasi-euclidean');
            dist_sum = dist1 + dist2;
            dist_sum = dist_sum - 1e8;
            if min(dist_sum(:)) < 1e8
                mask3 = dist_sum <= min(dist_sum(:))*1.0001;    % shortest path
                mask_final = zeros(x,y,z);
                mask_final(x_min:x_max, y_min:y_max, z_min:z_max) = mask3;
                
                g = addedge(g,ii,jj,min(dist_sum(:)));
                fore_candidate{size(g.Edges.EndNodes,1),1} = find(mask_final);      
            end
        end
    end
end
% save(fullfile(file_name, 'redetection_info_v2', ['t' num2str(tt) '_d' num2str(max_dist) '_before.mat']),'g','fore_candidate');
% plotShortestPath(foreground, g, fore_candidate, 20, ...
%      fullfile(file_name, 'redetection_info_v2', ['t' num2str(tt) '_d' num2str(max_dist) '_before']), conn);


%% min-span tree to remove redundant connections
% load(fullfile(file_name, 'redetection_info_v2', ['t' num2str(tt) '_d' num2str(max_dist) '_before.mat']));
num_node = height(g.Nodes);
G = zeros(num_node,num_node);
for kk = 1:size(g.Edges.EndNodes,1)
    ii = g.Edges.EndNodes(kk,1);
    jj = g.Edges.EndNodes(kk,2);
    G(ii,jj) = 1;
    G(jj,ii) = 1;
end
G = sparse(G);
[S,C] = graphconncomp(G);
g_new = graph;
for mm = 1:S
    if sum(C==mm) > 1
        node_list = find(C==mm);
        g_temp = graph;
        for kk = 1:size(g.Edges.EndNodes,1)
            ii = g.Edges.EndNodes(kk,1);
            jj = g.Edges.EndNodes(kk,2);
            if ismember(ii, node_list)
                % must use string, otherwise number of nodes will start from 1
                % and get wrong result
                g_temp = addedge(g_temp, num2str(ii), num2str(jj), g.Edges.Weight(kk));
            end
        end
        g_temp = minspantree(g_temp);
        for kk = 1:size(g_temp.Edges.EndNodes,1)
            if g_temp.Edges.Weight(kk) < priCvt_thre
                g_new = addedge(g_new, g_temp.Edges.EndNodes(kk,1), ...
                    g_temp.Edges.EndNodes(kk,2), g_temp.Edges.Weight(kk));
            end
        end
    end
end

fore_candidate_new = {};
% data_grad = imgradient3(imgaussfilt3(data,gaussian_sigma));
data_grad = imgradient3(data);
foreground_new = foreground;
for kk = 1:size(g_new.Edges.EndNodes,1)
%     kk
    ii = str2num(g_new.Edges.EndNodes{kk,1});
    jj = str2num(g_new.Edges.EndNodes{kk,2});
    [~, mm] = ismember([min(ii,jj) max(ii,jj)],g.Edges.EndNodes,'rows');   
%     fore_candidate_new{kk} = fore_candidate{mm};
    
    % min-cut segmentation
    [x_min, y_min, z_min, x_max, y_max, z_max] = getUnionRange(foreground, com, ii,jj);
    grad_temp = data_grad(x_min:x_max, y_min:y_max, z_min:z_max);
    foreground_temp = foreground(x_min:x_max, y_min:y_max, z_min:z_max);
    sMap = false(x,y,z);
%     sMap(com.PixelIdxList{ii}) = 1;
%     sMap(com.PixelIdxList{jj}) = 1;
    sMap(fore_candidate{mm}) = 1;
    tMap = true(x,y,z) - imdilate(sMap,strel("sphere",1));
    fMap = true(x,y,z);
    sMap = logical(sMap(x_min:x_max, y_min:y_max, z_min:z_max));
    fMap = logical(fMap(x_min:x_max, y_min:y_max, z_min:z_max));
    tMap = logical(tMap(x_min:x_max, y_min:y_max, z_min:z_max));
    
    % 6-conncection get much better segmentation. Don't use 26.
    [dat_in, src, sink] = graphCut_negHandle_mat(grad_temp, fMap, sMap, ...
        tMap, 6, [1 4], false);
    G = digraph(dat_in(:,1),dat_in(:,2),dat_in(:,3));
    if ~isempty(find(isnan(dat_in(:)), 1)) || isnan(sink)
        keyboard;
    end
    [~,~,cs,~] = maxflow(G, src, sink); % cs: fg, ct: bg
    cs = cs(cs<=numel(fMap));
    sMap_temp = zeros(size(sMap));
    sMap_temp(cs) = 1;
    sMap = false(x,y,z);
    sMap(x_min:x_max, y_min:y_max, z_min:z_max) = sMap_temp;
    sMap(com.PixelIdxList{ii}) = 0;
    sMap(com.PixelIdxList{jj}) = 0;
    fore_candidate_new{kk} = find(sMap);
    foreground_new(fore_candidate_new{kk}) = 1;
end
save(fullfile(file_name, 'redetection_info_v2', 'info', ['t' num2str(tt) '_d' num2str(max_dist) '.mat']),'g_new','fore_candidate_new');
plotShortestPath(foreground, g_new, fore_candidate_new, 2, ...
     fullfile(file_name, 'redetection_info_v2', 'connection', num2str(tt)), conn);
% remove tiny objbects
com = bwconncomp(foreground_new, conn);
for ii = 1:com.NumObjects
    if length(com.PixelIdxList{ii}) < 200
        foreground_new(com.PixelIdxList{ii}) = 0;
    end
end
tifwrite(uint8(foreground_new*255),  fullfile(file_name, 'redetection_info_v2', 'foreground', num2str(tt)));
end

function [x_min, y_min, z_min, x_max, y_max, z_max] = getUnionRange(foreground, com, ii,jj)
    [x,y,z] = size(foreground);
    [x_temp1, y_temp1, z_temp1] = ind2sub(size(foreground),com.PixelIdxList{ii});
    [x_temp2, y_temp2, z_temp2] = ind2sub(size(foreground),com.PixelIdxList{jj});
    x_min = min(min(x_temp1),min(x_temp2)); x_max = max(max(x_temp1),max(x_temp2)); 
    y_min = min(min(y_temp1),min(y_temp2)); y_max = max(max(y_temp1),max(y_temp2)); 
    z_min = min(min(z_temp1),min(z_temp2)); z_max = max(max(z_temp1),max(z_temp2)); 
    x_gap = x_max - x_min;
    y_gap = y_max - y_min;
    z_gap = z_max - z_min;
    x_min = max(1, round(x_min - 0.1*x_gap)); x_max = min(x, round(x_max + 0.1*x_gap));
    y_min = max(1, round(y_min - 0.1*y_gap)); y_max = min(y, round(y_max + 0.1*y_gap));
    z_min = max(1, round(z_min - 0.1*z_gap)); z_max = min(z, round(z_max + 0.1*z_gap));
end

function plotShortestPath(foreground, g, fore_candidate, threshold, path, conn)
%% plot the shortest path
[x,y,z] = size(foreground);
data_c = zeros(x,y,3,z);
data_r = zeros(x,y,z);
for ii = 1:numel(fore_candidate)
    data_temp = zeros(x,y,z);
    data_temp(fore_candidate{ii}) = 255; % - g.Edges.Weight(ii)*250/threshold;
%     data_temp(data_temp<0) = 0;
%     data_temp = imdilate(data_temp,conn);
    data_r = data_r + data_temp;
end
for zz = 1:z
    data_c(:,:,1,zz) = foreground(:,:,zz)*128 + data_r(:,:,zz);
    data_c(:,:,2,zz) = foreground(:,:,zz)*128;
    data_c(:,:,3,zz) = foreground(:,:,zz)*128;
end
tifwrite(uint8(data_c), path);
end

