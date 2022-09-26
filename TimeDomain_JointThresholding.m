clc;clear;close all;

%% get threshold result

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
mkdir([file_name '\foreground_iter1_jointThresholding']);
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
minsize = 3;
maxsize = 10000;
info_table = cell(t,1);
% OST
for tt = 1:2
    tt
    data = data_all(:,:,:,tt);
    z_map1 = nan(size(data));
    data_smoothed = imgaussfilt3(data,gaussian_sigma);

    x_direction = [1 -1 0 0 0 0];
    y_direction = [0 0 1 -1 0 0];
    z_direction = [0 0 0 0 1 -1];

    foreground = zeros(x,y,z);
    tic;

    data_nan = data;
    data_nan(foreground==1) = nan;
    score_map = zeros(x,y,z);
    z_map_temp = nan(x,y,z);
    components_ind = 0;
    components_cell = {};
    z_scores_list = [];
    threshold_list = [];
    tree_list = {};
    for thre = 1:25
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
                    components_ind = components_ind + 1;
                    components_cell{components_ind,1} = components.PixelIdxList{com_ind};
                    z_scores_list(components_ind,1) = z_score;
                    threshold_list(components_ind,1) = thre;
                    tree_temp = struct;
                    tree_temp.Father = [];
                    tree_temp.Children = [];
                    for rr = components_ind-1:-1:1
                        if all(ismember(components_cell{components_ind},components_cell{rr}))
                            tree_list{rr}.Children = [tree_list{rr}.Children components_ind];
                            tree_temp.Father = rr;
                            break;
                        end
                    end
                    tree_list{components_ind,1} = tree_temp;
                    
                    score_map(components.PixelIdxList{com_ind}) = max(z_score,score_map(components.PixelIdxList{com_ind}));
                    
                end
            end
        end
    end
    save([file_name '\foreground_iter1_jointThresholding\score_map_' num2str(tt)],'score_map');
%         components_cell{iter}(thre,1) = {components.PixelIdxList};
%         components_cell{iter}(thre,2) = {z_scores};
info_table{tt} = table(components_cell',z_scores_list',threshold_list',tree_list','VariableNames',{'Components','Zscore','Threshold','Tree'});
end
save([file_name '\foreground_iter1_jointThresholding\info_table'],'info_table');


%% use temporal results
clc;clear;close all;
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
method_name = 'foreground_iter1_jointThresholding5';
load([file_name '\data_c2']);
[x,y,z,t] = size(data2);
load([file_name '\' method_name '\info_table']);
%% initialization
t = 2; % temporary
time_length = t;
thre_length = 25;
thre_list = 10:10:250;
smoothness = 0.01;
com_num = zeros(time_length,1);
for tt = 1:time_length
    com_num(tt) = length(info_table{tt}.Zscore);
end
com_sum = [0 cumsum(com_num(1:end-1))];
s_ind = (thre_length-1)*sum(com_num) + 1;
t_ind = s_ind + 1;
%%
% build graph
% total node num: 24*com_num + 2
% each com node: (tt,ii) ->  1:24 + 24*(cum_sum(tt) + ii -1)
G = digraph();
% add sub graph
node_table = zeros( sum(com_num),2 );
for tt = 1:time_length
    for ii = 1:com_num(tt)
        s_list = [s_ind (1:thre_length-1) + 24*(com_sum(tt) + ii -1)];
        t_list = [(1:thre_length-1) + 24*(com_sum(tt) + ii -1) t_ind];
        weight = inf(1, thre_length);
        thre_temp = info_table{tt}.Threshold(ii);
        weight(thre_temp) = 1/info_table{tt}.Zscore(ii);
        jj = ii;
        while ~isempty(info_table{tt}.Tree{jj}.Father)
            jj = info_table{tt}.Tree{jj}.Father;
            thre_temp = info_table{tt}.Threshold(jj);
            weight(thre_temp) = 1/info_table{tt}.Zscore(jj);
        end
        G = addedge(G,s_list,t_list,weight);
        node_table(com_sum(tt) + ii,1) = tt;
        node_table(com_sum(tt) + ii,2) = ii;
    end
end
% add smoness edge
max_smo_num = sum(com_num(1:end-1).*com_num(2:end))*(thre_length-1);
count = 0;
% % before method 4
% for tt = 1:time_length-1
%     for ii = 1:com_num(tt)
%         ii
%         for jj = 1:com_num(tt+1)
%             list_ii = info_table{tt}.Components{ii};
%             list_jj = info_table{tt+1}.Components{jj};
%             if info_table{tt}.Threshold(ii) == info_table{tt+1}.Threshold(jj)  % only under the same threshold are connected
% %                 iou = length(intersect(list_ii,list_jj))/length(union(list_ii,list_jj));
%                 iou = length(intersect(list_ii,list_jj))/min(length(list_ii),length(list_jj)); % add more penalty on small components
%                 if iou > 0
%                     s_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt) + ii - 1);
%                     t_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt+1) + jj - 1);
%                     weight(count*(thre_length-1)+(1:thre_length-1)) = ones(1,24)*smoothness*iou;
%                     count = count + 1;
%                 end            
%             end
%         end
%     end
% end
% method 5
s_list = zeros(1,1000000);
t_list = zeros(1,1000000);
weight = zeros(1,1000000);
for tt = 1:time_length
    for ii = 1:com_num(tt)
        ii
        if tt < time_length
        for jj = 1:com_num(tt+1)
            list_ii = info_table{tt}.Components{ii};
            list_jj = info_table{tt+1}.Components{jj};
            if info_table{tt}.Threshold(ii) == info_table{tt+1}.Threshold(jj)  % only under the same threshold are connected
%                 iou = length(intersect(list_ii,list_jj))/length(union(list_ii,list_jj));
                iou = length(intersect(list_ii,list_jj))/min(length(list_ii),length(list_jj)); % add more penalty on small components
                if iou > 0
                    s_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt) + ii - 1);
                    t_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt+1) + jj - 1);
                    weight(count*(thre_length-1)+(1:thre_length-1)) = ones(1,24)*smoothness*iou;
                    count = count + 1;
                end            
            end
        end
        end
        % Connect parent and children
        for jj = 1:com_num(tt)
            list_ii = info_table{tt}.Components{ii};
            list_jj = info_table{tt}.Components{jj};
            if any(info_table{tt}.Tree{ii}.Father == jj) || any(info_table{tt}.Tree{jj}.Father == ii)  % only under the same threshold are connected
%                 iou = length(intersect(list_ii,list_jj))/min(length(list_ii),length(list_jj)); % add more penalty on small components
%                 if iou > 0
                    s_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt) + ii - 1);
                    t_list(count*(thre_length-1)+(1:thre_length-1)) = (1:thre_length-1) + 24*(com_sum(tt) + jj - 1);
                    weight(count*(thre_length-1)+(1:thre_length-1)) = ones(1,24)*smoothness;
                    count = count + 1;
%                 end            
            end
        end
    end
end
s_list = s_list(1:count*(thre_length-1));
t_list = t_list(1:count*(thre_length-1));
weight = weight(1:count*(thre_length-1));
G = addedge(G,s_list,t_list,weight);
G = addedge(G,t_list,s_list,weight);

%%
load([file_name '\' method_name '\graph_0.01.mat']);
% calculate max flow
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
                info_table{tt_temp}.Tree{ii_temp}.Cut = s_layer + 1;
            else
                info_table{tt_temp}.Tree{ii_temp}.Cut = t_layer;
            end
        end
    end
end
% get local maximum
local_max_flag = zeros(sum(com_num),1);
for tt = 1:time_length
    for ii = 1:com_num(tt)
        flag_temp = 0;
        if  info_table{tt}.Tree{ii}.Cut == info_table{tt}.Threshold(ii)
            flag_temp = 1;
            child_list = info_table{tt}.Tree{ii}.Children;
            while ~isempty(child_list)
                child_temp = child_list(1);
                child_list(1) = [];
                if info_table{tt}.Tree{child_temp}.Cut > info_table{tt}.Tree{ii}.Cut
                    flag_temp = 0;
                    break;
                end
                child_list = [child_list info_table{tt}.Tree{child_temp}.Children];
            end
        end
        local_max_flag(com_sum(tt) + ii) = flag_temp;
    end
end

%%
% output result
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);
for tt = 1:2 %time_length
    z_map1 = zeros(301,301,41);
    for ii = 1:com_num(tt)
        if local_max_flag(com_sum(tt) + ii) == 1
            z_map1(info_table{tt}.Components{ii}) = info_table{tt}.Zscore(ii);
        end
    end
    foreground = z_map1>2;
    % %% remove tiny objects
    components = bwconncomp(foreground,conn);
    foreground = zeros(size(foreground));
    for ii = 1:components.NumObjects
        if length(components.PixelIdxList{ii}) > 20
            foreground(components.PixelIdxList{ii}) = 1;
        end
    end
    tifwrite(uint8(foreground*255),[file_name '\' method_name '\temp_0.01_' num2str(tt)]);
end

%% validation
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);

% for ii = 1:s_cut
%     ss = s_cut(ii);
%     ind = 857;
%     if ss<= ind*24 && ss > (ind - 1)*24
%         if ismember(ss+1,t_cut)
%             mod(ss,24) + 1
%         end
%     end
% end
% 
% for ii = 1:4441
%     if ismember(1133742,info_table{tt}.Components{ii})
%         fprintf('%d:',ii);
% %         if local_max_flag(ii) == 1
% %             fprintf('yes\n');
% %         else
% %             fprintf('no\n');
% %         end
%         fprintf('%d %d\n', info_table{tt}.Threshold(ii), info_table{tt}.Tree{ii}.Cut);
%     end
% end
for tt = 1:1 %time_length
    z_map1 = zeros(301,301,41);
    for ii = 1:com_num(tt)
        if local_max_flag(com_sum(tt) + ii) == 1
            z_map1(info_table{tt}.Components{ii}) = info_table{tt}.Zscore(ii);
        end
    end
    z_map1 = z_map1>0;
%     z_map1 = imregionalmax(z_map1,26);
%     components = bwconncomp(z_map1,conn);
%     components.NumObjects
%     for ii = 1:com_num(tt)
%         if local_max_flag(com_sum(tt) + ii) == 1
%             flag = 0;
%             for jj = 1:components.NumObjects
% %                 if length(info_table{tt}.Components{ii}) == length(components.PixelIdxList{jj}) ...
% %                         && all(info_table{tt}.Components{ii} == components.PixelIdxList{jj})
%                 if all(ismember(info_table{tt}.Components{ii},components.PixelIdxList{jj}))
%                     flag = 1;
%                     break
%                 end
%             end
%             if ~flag
%                 ii
%             end
% %             if ismember(3563554,info_table{tt}.Components{ii})
% %                 ii
% %             end
%         end
%     end

%     for jj = 1:components.NumObjects
%         if any(ismember(3563554,components.PixelIdxList{jj}))
%             jj
%             break
%         end
%     end



    z_map2 = zeros(301,301,41);
    for ii = 1:com_num(tt)
        z_map2(info_table{tt}.Components{ii}) = max(z_map2(info_table{tt}.Components{ii}), info_table{tt}.Zscore(ii));
    end
    z_map2 = imregionalmax(z_map2,conn);
    components = bwconncomp(z_map2,conn);
    components.NumObjects

    
    im = zeros(301,301,3,41);
    for zz = 1:41
        im(:,:,1,zz) = z_map2(:,:,zz);
        im(:,:,2,zz) = z_map1(:,:,zz);
        im(:,:,3,zz) = z_map1(:,:,zz)&z_map2(:,:,zz);
    end
    tifwrite(uint8(im*255),[file_name '\' method_name '\temp' num2str(tt)]);
end

%% temp result check
% x_loc = 151; y_loc = 213; z_loc = 38;  1/0
x_loc = 25; y_loc =54; z_loc = 28;  % 0/1
loc = sub2ind([301 301 41],x_loc,y_loc,z_loc);
loc = 2280693;
% load([file_name '\foreground_iter1_jointThresholding\info_table']);
% load([file_name '\foreground_iter1_jointThresholding\graph_0.01.mat']);
% calculate max flow
[~,~,s_cut,t_cut] = maxflow(G,s_ind,t_ind);
EndNodes = G.Edges.EndNodes;
% get cut threshold for every component
temp = [];
for tt = 1:1 %time_length
    for ii = 1:com_num(tt)
        if ismember(loc,info_table{tt}.Components{ii})
            temp = [temp ii];
        end
    end
end
