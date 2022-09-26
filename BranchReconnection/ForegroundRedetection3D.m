%% Initialization
clc;clear;close all;
dbstop if error
% it will re-detect foreground based on the order-statistic thresholding
% result, connect broken branches and remove connected tips

% Algorithm:
% Input: foreground (from order-stattistic thresholding or sth else)
% Parameter: distance d
% 1. For every two component whose distance smaller than d:
%       1.1 Find the shorest path
%       1.2 Calculate the z-score based on hypothest testing
% 2. Use the minimum spanning tree to pune abudnant results

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';
mkdir([file_name 'redetection_info']);
data_all = load([file_name 'stabilization_c2\stabilization_data.mat']);
data_all = data_all.data;
[x,y,z,t] = size(data_all);

%% tt = 20;
time1 = 0;
time2 = 0;
time3 = 0;
tic;
for tt = 1:t
    tt
ind = num2str(1000+tt);
ind = ind(2:4);
foreground = load([file_name 'foreground_c2\' ind '.mat']);
foreground = foreground.foreground;

gaussian_sigma = 0.4;
d = 20;
z_map = zeros(x,y,z);
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);

com_before = bwconncomp(foreground,conn);
data = data_all(:,:,:,tt);
data = imgaussfilt3(data,gaussian_sigma);
cost_map = 1./(data+1);
cost_map(logical(foreground)) = 1;

% stage 1
% % pick index    % 85 115 / 120 126
% for ii = 1:com_before.NumObjects
%     com_ind = sub2ind(size(data),152,161,12);
%     if ismember(com_ind,com_before.PixelIdxList{ii})
%         ii
%     end
% end
% candidate
g = graph;
fore_candidate = cell(1,1);
for ii = 1:com_before.NumObjects
%     ii
    mask_temp = zeros(x,y,z);
    mask_temp(com_before.PixelIdxList{ii}) = 1;
    dist_temp = bwdist(mask_temp);
    
    for jj = ii+1:com_before.NumObjects
        if min(dist_temp(com_before.PixelIdxList{jj})) < d && jj~=ii
            
            time_start = toc;
            % find shortest path
            % mask1&2: two foreground %mask3: candidate %mask4: neighbor in foreground
            mask1 = false(x,y,z);
            mask1(com_before.PixelIdxList{ii}) = 1;
            dist1 = graydist(cost_map, mask1, 'quasi-euclidean');
            mask2 = false(x,y,z);
            mask2(com_before.PixelIdxList{jj}) = 1;
            dist2 = graydist(cost_map, mask2, 'quasi-euclidean');
            dist_sum = dist1 + dist2;
            dist_sum = dist_sum - 1;
            dist_sum(logical(foreground)) = inf;
            mask3 = dist_sum <= min(dist_sum(:))*1.001;    % shortest path
%             dist3 = bwdist(mask3, 'quasi-euclidean');
            time1 = time1 + toc - time_start;
            
            time_start = toc;
            % caclulate z-score and get min span tree
            mask3 = imdilate(mask3,conn) & (~mask1&~mask2);
%             thre = median(data(mask3));
%             dist3 = bwdistgeodesic(data>thre, mask3);
%             mask3 = (dist3 < inf) | mask3;
            mask3 = regionGrowing(data, logical(foreground), mask3, conn) | mask3;
            mask4 = imdilate(mask3,ones(3,3,3));
            mask4 = imdilate(mask4,ones(3,3,3)) & (~logical(foreground)&~mask3&~mask4);
%             mask3 = findMinCut(data, foreground, mask1, mask2, mask3);
%             z_score = order_stat_test(sum(mask4(:)), sum(mask3(:)), mean(data(mask4)) - mean(data(mask3)));
            time2 = time2 + toc - time_start;
            
            time_start = toc;
            background = data(mask4);
            if length(background) > 4*sum(mask3(:))
                background = uniformSampling(background, 4*sum(mask3(:)));
            end
           
                
            [z_score, mu, sigma] = order_stat_test_boyu(data(mask3),background);


            if z_score > 2
                 g = addedge(g,ii,jj,z_score);
                 fore_candidate{size(g.Edges.EndNodes,1),1} = find(mask3);
            end
            time3 = time3 + toc - time_start;
       end
    end
%     toc
%     g.Edges.EndNodes
end
time1
time2
time3
toc
% save(['redetection_graph_t' num2str(tt) '_d' num2str(d) '_boyu_z2_40_new.mat'],'g','fore_candidate')
save([file_name 'redetection_info\t' num2str(tt) '_d' num2str(d) '.mat'],'g','fore_candidate')
end

%% Minspance tree & plot
tic;
% load redetection_graph_t20_d20_boyu_z2_40_new.mat
mkdir([file_name '\foreground_c2_reconnect']);
for tt = 33:35
ind = num2str(1000+tt);
ind = ind(2:4);
load([file_name 'redetection_info\t' num2str(tt) '_d20.mat']);
foreground = load([file_name 'foreground_c2\' ind '.mat']);
foreground = foreground.foreground;
gaussian_sigma = 0.4;
data = data_all(:,:,:,tt);
data = imgaussfilt3(data,gaussian_sigma);
cost_map = 1./(data+1);
cost_map(logical(foreground)) = 1;



% get the min span tree
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
%                 g_temp = addedge(g_temp, num2str(ii), num2str(jj), g.Edges.Weight(kk));
                g_temp = addedge(g_temp, num2str(ii), num2str(jj), sum(cost_map(fore_candidate{kk})));
            end
        end
        g_temp = minspantree(g_temp);
        g_new = addedge(g_new, g_temp.Edges);
        a = 1;
    end
end
% fore_temp = foreground*128;
fore_temp = zeros(x,y,3,z);
fore_new = foreground;
for zz = 1:z
    fore_temp(:,:,1,zz) = foreground(:,:,zz)*128;
    fore_temp(:,:,2,zz) = foreground(:,:,zz)*128;
    fore_temp(:,:,3,zz) = foreground(:,:,zz)*128;
end
for kk = 1:size(g_new.Edges.EndNodes,1)
    kk
    ii = str2num(g_new.Edges.EndNodes{kk,1});
    jj = str2num(g_new.Edges.EndNodes{kk,2});
    
    % find shortest path
    mask3 = false(x,y,z);
    [~, mm] = ismember([min(ii,jj) max(ii,jj)],g.Edges.EndNodes,'rows');
    mask3(fore_candidate{mm}) = 1;
    fore_new = fore_new + mask3;
    color = rand(3,1)/2+0.5;
    for zz = 1:z
        fore_temp(:,:,1,zz) = mask3(:,:,zz)*255*color(1) + (1-mask3(:,:,zz)).*fore_temp(:,:,1,zz);
        fore_temp(:,:,2,zz) = mask3(:,:,zz)*255*color(2) + (1-mask3(:,:,zz)).*fore_temp(:,:,2,zz);
        fore_temp(:,:,3,zz) = mask3(:,:,zz)*255*color(3) + (1-mask3(:,:,zz)).*fore_temp(:,:,3,zz);
    end
    
end
foreground = fore_new;
% save('020','foreground');
% ind = num2str(1000+tt); 
% ind = [file_name '\reconnect_temp\' ind(2:4)];
% tifwrite(uint8(fore_temp),ind);
% tifwrite(uint8(fore_temp),'temp1');
% tifwrite(uint8(fore_new*255),'temp2');

ind = num2str(1000+tt); 
ind = [file_name '\foreground_c2_reconnect\' ind(2:4)];
save(ind,'foreground');
tifwrite(uint8(foreground*255),ind);
    
end
toc

%% plot

load redetection_graph_t20_d20_boyu_z2_old.mat
fore_temp = zeros(x,y,3,z);
for zz = 1:z
    fore_temp(:,:,1,zz) = foreground(:,:,zz)*128;
    fore_temp(:,:,2,zz) = foreground(:,:,zz)*128;
    fore_temp(:,:,3,zz) = foreground(:,:,zz)*128;
end
for kk = 1:size(g.Edges.EndNodes,1)
    kk
    ii = g.Edges.EndNodes(kk,1);
    jj = g.Edges.EndNodes(kk,2);
    
    % find shortest path
    mask3 = false(x,y,z);
    mask3(fore_candidate{kk}) = 1;
    color = rand(3,1)/2+0.5;
    for zz = 1:z
        fore_temp(:,:,1,zz) = mask3(:,:,zz)*255*color(1) + (1-mask3(:,:,zz)).*fore_temp(:,:,1,zz);
        fore_temp(:,:,2,zz) = mask3(:,:,zz)*255*color(2) + (1-mask3(:,:,zz)).*fore_temp(:,:,2,zz);
        fore_temp(:,:,3,zz) = mask3(:,:,zz)*255*color(3) + (1-mask3(:,:,zz)).*fore_temp(:,:,3,zz);
    end
    
end
tifwrite(uint8(fore_temp),'temp2')



function result = uniformSampling(data, m)
    data = sort(data);
    n = length(data);
    if m > n
        error('Sampling size is larger than data size.');
    end
    loc = [1:m]*(n+1)/(m+1);
    result = interp1(data,loc);
end

function mask3 = regionGrowing(data, foreground, mask3, conn)
    % growing region until no satisfactory
%     thre = median(data(mask3));
    thre = prctile(data(mask3), 40);
    mask_bound = imdilate(mask3,ones(3,3,3));
%     mask_bound = imdilate(mask_bound,ones(3,3,3));
    mask_old = false(size(mask3));
    while any(xor(mask3, mask_old),'all')
        mask_old = mask3;
        mask3 = imdilate(mask3,conn);
        mask3(data<thre) = 0;
        mask3(foreground) = 0;
    end
    mask3 = mask3.*mask_bound;
end

function mask3 = findMinCut(data, foreground, mask1, mask2, mask3)
    x_direction = [ 0 -1 1  0 0 0];
    y_direction = [ 0  0 0 -1 1 0];
    z_direction = [-1  0 0  0 0 1];

    source_list = find(mask3);
%     mask_fore = mask1 + mask2;
%     dist_fore = bwdist(mask_fore, 'quasi-euclidean');
    dist3 = bwdist(mask3, 'quasi-euclidean');
    mask_candi = (dist3<5) & (~foreground&~mask3);
    mask_sink = (dist3<7) & (~foreground&~mask3&~mask_candi);
    candi_list = find(mask_candi);
    sink_list = find(mask_sink);
    g = graph;
    for ii = 1:length(source_list)
        g = addedge(g, 's', num2str(source_list(ii)), inf);
        [x_ind, y_ind, z_ind] = ind2sub(size(data), source_list(ii));
        for dd = 1:6
            try
                neighbor_ind = sub2ind(size(data),x_ind+x_direction(dd), y_ind+y_direction(dd), z_ind+z_direction(dd));      
                if ismember(neighbor_ind, candi_list)
                    g = addedge(g, num2str(source_list(ii)), num2str(neighbor_ind), 1/(abs(data(source_list(ii)) - data(neighbor_ind))+1));
                end
            catch ME
                if ~strcmp(ME.identifier,  'MATLAB:sub2ind:IndexOutOfRange')
                    error(ME.message);
                end
            end
        end
    end
    for ii = 1:length(sink_list)
        g = addedge(g, 't', num2str(sink_list(ii)), inf);
        [x_ind, y_ind, z_ind] = ind2sub(size(data), sink_list(ii));
        for dd = 1:6
            try
                neighbor_ind = sub2ind(size(data),x_ind+x_direction(dd), y_ind+y_direction(dd), z_ind+z_direction(dd));
                if ismember(neighbor_ind, candi_list)
                    g = addedge(g, num2str(sink_list(ii)), num2str(neighbor_ind), 1/(abs(data(sink_list(ii)) - data(neighbor_ind))+1));
                end
            catch ME
                if ~strcmp(ME.identifier,  'MATLAB:sub2ind:IndexOutOfRange')
                    error(ME.message);
                end
            end
        end
    end
    for ii = 1:length(candi_list)
        [x_ind, y_ind, z_ind] = ind2sub(size(data), candi_list(ii));
        for dd = 1:6
            try
                neighbor_ind = sub2ind(size(data),x_ind+x_direction(dd), y_ind+y_direction(dd), z_ind+z_direction(dd));
                if ismember(neighbor_ind, candi_list)
                    g = addedge(g, num2str(candi_list(ii)), num2str(neighbor_ind), 1/(abs(data(candi_list(ii)) - data(neighbor_ind))+1));
                end
            catch ME
                if ~strcmp(ME.identifier,  'MATLAB:sub2ind:IndexOutOfRange')
                    error(ME.message);
                end
            end
        end
    end
    [mf,GF,cs,ct] = maxflow(g, 's', 't');
    mask4 = mask3;
    for ii = 1:length(cs)
        if ~strcmp(cs{ii}, 's') && ismember(str2num(cs{ii}), candi_list)
            mask3(str2num(cs{ii})) = 1;
        end
    end
end

function z_score = order_stat_test(M,N,L)
    pdfz = normpdf(norminv(N/(M+N),0,1));
    cdfz = normcdf(norminv(N/(M+N),0,1));
    mu = pdfz/(1-cdfz) + pdfz/cdfz;
    sigma = sqrt((1+norminv(N/(M+N),0,1)*pdfz/(1-cdfz)-(pdfz/(1-cdfz))^2)/M + (1-norminv(N/(M+N),0,1)*pdfz/cdfz-(pdfz/cdfz)^2)/N);
    z_score = (L/sqrt(767) - mu)/sigma;  %variance for certian data
end

function [z_score, mu, sigma] = order_stat_test_boyu(largeGroup, smallGroup)
    [mu, sigma] = ksegments_orderstatistics_v2(largeGroup, smallGroup);
    L = mean(largeGroup) - mean(smallGroup);
    z_score = (L/sqrt(767) - mu)/sigma;  %variance for certian data
end

