clc;clear;close all;

%% get threshold result

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
method_name = 'foreground_iter1_jointThresholding2';
mkdir([file_name '\' method_name]);
load([file_name '\data_c2']);
load([file_name '\stabilizeFunction_c2']);
data_all = double(data2);
[x,y,z,t] = size(data_all);
std = sqrt(max(variance)/3);

data_all = interp1(0:255,stabilizeFunction,data_all);

% get data
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

% 6 neighbor connecivty
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);
gaussian_sigma = 0.8;
z_thre = 2;
minsize = 5;
maxsize = 30000;
info_table = cell(t,1);
% OST
for tt = 1:t
    tt
    data = data_all(:,:,:,tt);
    z_map1 = nan(size(data));
    data_smoothed = imgaussfilt3(data,gaussian_sigma);
    score_com = cell(25,1);
    
    x_direction = [1 -1 0 0 0 0];
    y_direction = [0 0 1 -1 0 0];
    z_direction = [0 0 0 0 1 -1];

    foreground = zeros(x,y,z);
    tic;

    data_nan = data;
    data_nan(foreground==1) = nan;
    score_map = zeros(x,y,z,25);
    components_ind = 0;
    components_cell = {};
    z_scores_list = [];
    threshold_list = [];
    tree_list = {};
    for thre = 1:25
        score_map_temp = zeros(x,y,z);
        binaryData = zeros(x,y,z);
        binaryData(data_smoothed<= thre*10) = 0;
        binaryData(data_smoothed>  thre*10) = 1;
        binaryData(foreground == 1) = 0;
        components = bwconncomp(binaryData,conn);
        components_len = length(components.PixelIdxList);
        for com_ind = 1:components_len
            com_element = components.PixelIdxList{com_ind};
            if length(com_element)>= minsize && length(com_element) <= maxsize
                [x_position,y_position,z_position] = ind2sub([x,y,z],com_element);
                com_neighbor = [];
                for dir_ind = 1:length(x_direction)
                    x_neighbor = min(max(x_position + x_direction(dir_ind),1),x); 
                    y_neighbor = min(max(y_position + y_direction(dir_ind),1),y); 
                    z_neighbor = min(max(z_position + z_direction(dir_ind),1),z); 
                    pix1 = sub2ind([x,y,z],x_neighbor,y_neighbor,z_neighbor);
                    pix1 = pix1(~ismember(pix1,com_element));
                    pix1 = pix1(~ismember(pix1,com_neighbor));
                    com_neighbor = [com_neighbor;pix1];
                end
                com_neighbor = unique(com_neighbor);
                com_neighbor(isnan(data_nan(com_neighbor))) = [];
                
                L = mean(data_nan(com_element)) - mean(data_nan(com_neighbor));
                largeGroup = data_nan(com_element);
                smallGroup = data_nan(com_neighbor);
                [mu, sigma] = ordStatApproxKsec(largeGroup, smallGroup);
                z_score = (L/std - mu)/sigma;
                if z_score > 0
%                     components_ind = components_ind + 1;
%                     components_cell{components_ind,1} = components.PixelIdxList{com_ind};
%                     z_scores_list(components_ind,1) = z_score;
%                     threshold_list(components_ind,1) = thre;
%                     tree_temp = struct;
%                     tree_temp.Father = [];
%                     tree_temp.Children = [];
%                     for rr = components_ind-1:-1:1
%                         if all(ismember(components_cell{components_ind},components_cell{rr}))
%                             tree_list{rr}.Children = [tree_list{rr}.Children components_ind];
%                             tree_temp.Father = rr;
%                             break;
%                         end
%                     end
%                     tree_list{components_ind,1} = tree_temp;
                    
                    score_map_temp(components.PixelIdxList{com_ind}) = z_score;
                    score_com{thre}{end+1} = components.PixelIdxList{com_ind};
                end
            end
        end
        score_map(:,:,:,thre) = score_map_temp;
    end
    save([file_name '\' method_name '\score_map_' num2str(tt)],'score_map','score_com');

% info_table{tt} = table(components_cell',z_scores_list',threshold_list',tree_list','VariableNames',{'Components','Zscore','Threshold','Tree'});
end
% save([file_name '\' method_name '\info_table'],'info_table');
%%
clear;
% get local maximum and component tree
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
method_name = 'foreground_iter1_jointThresholding2';

t = 33; % temporary
time_length = t;
thre_length = 25;
thre_list = 10:10:250;
thre_dist = 4;
smoothness = 0.1;
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);
%%
info_table = cell(t,1);
for tt = 1:t
tt
load([file_name '\' method_name '\score_map_' num2str(tt+1)]);
[score_max, score_max_ind] = max(score_map,[],4);
score_max_ind(score_max==0) = 0;

% get relationship
foreground = imregionalmax(score_max,conn);
foreground = foreground > 0;
fore_com = bwconncomp(foreground,conn);

thre_current_list = zeros(fore_com.NumObjects,1);
thre_candidate_list = cell(fore_com.NumObjects,1);
z_score_list = cell(fore_com.NumObjects,1);
pixel_list = cell(fore_com.NumObjects,1);
for ii = 1:fore_com.NumObjects
    % check z_score unique
    if length(unique(score_max(fore_com.PixelIdxList{ii}))) > 1
        error('Not unique zscore');
    end
    if length(unique(score_max_ind(fore_com.PixelIdxList{ii}))) > 1
        error('Not unique region');
    end
    
    thre_current = unique(score_max_ind(fore_com.PixelIdxList{ii}));
    thre_candidate = max(1,thre_current-thre_dist):min(thre_length,thre_current+thre_dist);
    z_score = zeros(length(thre_candidate),1);
    pixel = cell(length(thre_candidate),1);
    for thre = max(1,thre_current-thre_dist):min(thre_length,thre_current+thre_dist)
        jj = thre - max(1,thre_current-thre_dist) + 1;
        if thre < thre_current
            score_map_temp = score_map(:,:,:,thre);
            for kk = 1:length(score_com{thre})
%                 if any(ismember(fore_com.PixelIdxList{ii}, score_com{thre}{kk}))
%                     all(ismember(fore_com.PixelIdxList{ii}, score_com{thre}{kk}))
% %                     error('Parent dont contian');
%                 end
                if all(ismember(fore_com.PixelIdxList{ii}, score_com{thre}{kk}))
                    if length(unique(score_map_temp(score_com{thre}{kk}))) > 1
                        error('multi scores');
                    end
                    z_score(jj) = unique(score_map_temp(score_com{thre}{kk}));
                    pixel{jj} = score_com{thre}{kk};
                    break;
                end
            end
        end
        if thre == thre_current
            z_score(jj) = unique(score_max(fore_com.PixelIdxList{ii}));
            pixel{jj} = fore_com.PixelIdxList{ii};
        end
        if thre > thre_current
            score_map_temp = score_map(:,:,:,thre);
            for kk = 1:length(score_com{thre})
                if all(ismember(score_com{thre}{kk}, fore_com.PixelIdxList{ii}))
                    if length(unique(score_map_temp(score_com{thre}{kk}))) > 1
                        error('multi scores');
                    end
                    z_score(jj) = max(z_score(jj), unique(score_map_temp(score_com{thre}{kk})));
                    pixel{jj} = [pixel{jj}; score_com{thre}{kk}];
                end
            end
        end
        
    end
    thre_current_list(ii) = thre_current;
    thre_candidate_list{ii} = thre_candidate;
    z_score_list{ii} = z_score;
    pixel_list{ii} = pixel;
end
result = zeros(size(thre_current_list));
info_table{tt} = table(pixel_list,z_score_list,thre_current_list,thre_candidate_list,result,'VariableNames',{'Components','Zscore','Threshold','Candidate','Result'});
end

%%  build graph
com_num = cellfun(@(x)size(x,1), info_table);
com_sum = [0; cumsum(com_num(1:end-1))];
s_ind = (thre_length-1)*sum(com_num) + 1;
t_ind = s_ind + 1;
node_table = zeros(sum(com_num),2);
G = digraph();
for tt = 1:time_length
    for ii = 1:com_num(tt)
        s_list = [s_ind (1:24) + 24*(com_sum(tt) + ii -1)];
        t_list = [(1:24) + 24*(com_sum(tt) + ii -1) t_ind];
        weight = inf(1,25);
        for jj = 1:length(info_table{tt}.Candidate{ii})
            weight(info_table{tt}.Candidate {ii}(jj)) = 1/info_table{tt}.Zscore{ii}(jj);
        end
        G = addedge(G,s_list,t_list,weight);
        node_table(com_sum(tt) + ii,1) = tt;
        node_table(com_sum(tt) + ii,2) = ii;
    end
end
% add smooth edge
s_list = zeros(1,1000000);
t_list = zeros(1,1000000);
weight = zeros(1,1000000);
count = 0;
for tt = 1:time_length-1
    tt
    for ii = 1:com_num(tt)
%         ii
        for jj = 1:com_num(tt+1)
            % find max iou
            list_ii = info_table{tt}.Components{ii};
            list_jj = info_table{tt+1}.Components{jj};
            iou = 0;
            for kk = 1:length(list_ii)
                for ll = 1:length(list_jj)
%                     iou = max(iou, length(intersect(list_ii{kk},list_jj{ll}))/length(union(list_ii{kk},list_jj{ll})));
                    iou = max(iou, length(intersect(list_ii{kk},list_jj{ll}))/min(length(list_ii{kk}),length(list_jj{ll})));
                end
            end
 
            if iou > 0
                s_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt) + ii - 1);
                t_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt+1) + jj - 1);
                weight(count*(thre_length-1)+(1:thre_length-1)) = ones(1,24)*smoothness*iou;
                count = count + 1;
            end            
        end
    end
end
s_list = s_list(1:count*(thre_length-1));
t_list = t_list(1:count*(thre_length-1));
weight = weight(1:count*(thre_length-1));
G = addedge(G,s_list,t_list,weight);
G = addedge(G,t_list,s_list,weight);

%% calculate max flow
[~,~,s_cut,t_cut] = maxflow(G,s_ind,t_ind);
count = 0;
EndNodes = G.Edges.EndNodes;
% get cut threshold for every component
for ii = 1:length(EndNodes)
    % _node: node index in graph
    % _layer: node index in subgraph
    % com_temp: component index
    s_node = EndNodes(ii,1);
    t_node = EndNodes(ii,2);
    if s_node == s_ind || t_node == t_ind || t_node == s_node + 1
        if ismember(s_node,s_cut) && ismember(t_node,t_cut)
            count = count + 1;
            
            if s_node ~= s_ind
                s_layer = mod(s_node,thre_length-1);
                com_temp = (s_node-s_layer)/(thre_length-1) + 1;
            else
                t_layer = mod(t_node,thre_length-1);
                com_temp = (t_node-t_layer)/(thre_length-1) + 1;
            end

            tt_temp = node_table(com_temp,1);
            ii_temp = node_table(com_temp,2);
            if s_node ~= s_ind
                info_table{tt_temp}.Result(ii_temp) = s_layer + 1;
            else
                info_table{tt_temp}.Result(ii_temp) = t_layer;
            end
        end
    end
end
%% comparison
for tt = 1:time_length

conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);
          
load([file_name '\' method_name '\score_map_' num2str(tt+1)]);
[score_max, score_max_ind] = max(score_map,[],4);
score_max_ind(score_max==0) = 0;

% get relationship
foreground = imregionalmax(score_max,conn);
% foreground = foreground > 0;

result = zeros(size(foreground));

for ii = 1:com_num(tt)
    result_temp = info_table{tt}.Components{ii}{info_table{tt}.Result(ii) - info_table{tt}.Candidate{ii}(1) + 1};
    result(result_temp) = 1;
    %     result() = 1;% ...
%         info_table{tt}.Zscore{ii}(info_table{tt}.Result(ii) - info_table{tt}.Candidate{ii}(1) + 1);
end
% score_max(~foreground) = 0;
% sum(abs(result-foreground),'all')
tifwrite(uint8(result*255),[file_name '\' method_name '\0.1_4\' num2str(tt+1)]);
end
% tifwrite(uint8(im*255),[file_name '\' method_name '\0.01']);
%% try to find all foreground
for tt = 1:time_length

conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);
          
load([file_name '\' method_name '\score_map_' num2str(tt+1)]);
result = zeros(size(score_map(:,:,:,1)));
for ii = 1:com_num(tt)
    result_temp = info_table{tt}.Components{ii}{info_table{tt}.Result(ii) - info_table{tt}.Candidate{ii}(1) + 1};
%     result(result_temp) = 1;
    layer_temp = info_table{tt}.Result(ii);
    if layer_temp == 1
        result(result_temp) = 1;
    else
        com_parent = bwconncomp(score_map(:,:,:,layer_temp-1),conn);
        find_parent_flag = 0;
        for jj = 1:com_parent.NumObjects
            if all(ismember(result_temp,com_parent.PixelIdxList{jj}))
                find_parent_flag = 1;
                break;
            end
        end
        if ~find_parent_flag
            result(result_temp) = 1;
        else
            score_map_child = score_map(:,:,:,layer_temp);
            com_child = bwconncomp(score_map(:,:,:,layer_temp),conn);
            for kk = 1:com_child.NumObjects
                if all(ismember(com_child.PixelIdxList{kk},com_parent.PixelIdxList{jj})) && unique(score_map_child(com_child.PixelIdxList{kk})) > 10
                    result(com_child.PixelIdxList{kk}) = 1;
                end
            end
        end
    end
end
% score_max(~foreground) = 0;
% sum(abs(result-foreground),'all')
tifwrite(uint8(result*255),[file_name '\' method_name '\0.1_4_union\' num2str(tt+1)]);
end
%% temp result check
x_loc = 283; y_loc = 170; z_loc = 2;
loc1 = sub2ind([301 301 41],x_loc,y_loc,z_loc);
x_loc = 287; y_loc = 131; z_loc = 14;
loc2 = sub2ind([301 301 41],x_loc,y_loc,z_loc);

tt = 3;
load([file_name '\' method_name '\score_map_' num2str(tt+1)]);
for thre = 1:25
    for ii = 1:length(score_com{thre})
        if all(ismember([loc1 loc2], score_com{thre}{ii}))
            fprintf('%d %d %d\n',thre,ii,length(score_com{thre}{ii}));
        end
    end
end

% [score_max, score_max_ind] = max(score_map,[],4);
% foreground = imregionalmax(score_max,conn);
% fore_com = bwconncomp(foreground,conn);

for ii = 1:length(info_table{tt}.Components)
    for jj = 1:length(info_table{tt}.Components{ii})
        if ismember(loc2, info_table{tt}.Components{ii}{jj})
            fprintf('%d %d %d\n',ii,jj,length(info_table{tt}.Components{ii}{jj}));
        end
    end
end
