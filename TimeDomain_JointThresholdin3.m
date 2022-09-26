clc;clear;
% get local maximum and component tree
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
method_name = 'foreground_iter1_jointThresholding3';
          
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
gaussian_sigma = 0.8;


t = 33; % temporary
time_length = t;
thre_length = 25;
thre_list = 10:10:250;
thre_dist = 2;
smoothness = 0.05;
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);

info_table = cell(t,1);
for tt = 1:time_length
% tt = 2;
tt
data = data_all(:,:,:,tt+1);
data_smoothed = imgaussfilt3(data,gaussian_sigma);
load([file_name '\foreground_iter1_jointThresholding2\score_map_' num2str(tt+1)]);
[score_max, score_max_ind] = max(score_map,[],4);
score_max_ind(score_max==0) = 0;
foreground = imregionalmax(score_max,conn);
fore_com = bwconncomp(foreground,conn);

thre_current_list = zeros(fore_com.NumObjects,1);
thre_candidate_list = cell(fore_com.NumObjects,1);
z_score_list = cell(fore_com.NumObjects,1);
pixel_list = cell(fore_com.NumObjects,1);

% region growing
x_direction = [1 -1 0 0 0 0];
y_direction = [0 0 1 -1 0 0];
z_direction = [0 0 0 0 1 -1];
[x,y,z] = size(score_max);
z_thre = 5;
fore_new = zeros(x,y,z);
for ii = 1:fore_com.NumObjects
    thre_current = unique(score_max_ind(fore_com.PixelIdxList{ii}));
    thre_candidate = max(1,thre_current-thre_dist):min(thre_length,thre_current+thre_dist);
    z_score = zeros(length(thre_candidate),1);
    pixel = cell(length(thre_candidate),1);
    
    com_neighbor = fore_com.PixelIdxList{ii};
    com_element = com_neighbor;
    while ~isempty(com_neighbor)
        [x_position,y_position,z_position] = ind2sub([x,y,z],com_neighbor);
        com_temp = [];
        for dir_ind = 1:length(x_direction)
            x_neighbor = min(max(x_position + x_direction(dir_ind),1),x); 
            y_neighbor = min(max(y_position + y_direction(dir_ind),1),y); 
            z_neighbor = min(max(z_position + z_direction(dir_ind),1),z); 
            pix1 = sub2ind([x,y,z],x_neighbor,y_neighbor,z_neighbor); % find all neighbor
            pix1(score_max(pix1) > score_max(com_neighbor)) = []; % remove larger pixel
            pix1 = pix1(~ismember(pix1,com_element)); % remove inside pixel
            pix1 = pix1(~ismember(pix1,com_temp)); % remove already selected pixel
            pix1(score_max(pix1) < z_thre) = [];      % remove insignificant pixel
            com_temp = [com_temp;pix1];
        end
        com_neighbor = com_temp;
        com_element = [com_element;com_neighbor];
    end
    fore_new(com_element) = 1;
    
    thre_current = unique(score_max_ind(fore_com.PixelIdxList{ii}));
    thre_candidate = max(1,thre_current-thre_dist):min(thre_length,thre_current+thre_dist);
    z_score = zeros(length(thre_candidate),1);
    pixel = cell(length(thre_candidate),1);
    for thre = max(1,thre_current-thre_dist):min(thre_length,thre_current+thre_dist)
        jj = thre - max(1,thre_current-thre_dist) + 1;
        pixel{jj} = com_element(data_smoothed(com_element) > thre*10);
        score_map_temp = score_map(:,:,:,thre);
        if ~isempty(pixel{jj})
            z_score(jj) = max(score_map_temp(pixel{jj}));
        end
    end
    
    thre_current_list(ii) = thre_current;
    thre_candidate_list{ii} = thre_candidate;
    z_score_list{ii} = z_score;
    pixel_list{ii} = pixel;
end
% tifwrite(uint8(fore_new*255),[file_name '\' method_name '\test\' num2str(tt+1)]);
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
mkdir([file_name '\' method_name '\' num2str(smoothness) '_' num2str(z_thre)]);
for tt = 1:time_length

% conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
%               [0 1 0; 1 1 1; 0 1 0],...
%               [0 0 0; 0 1 0; 0 0 0]);
%           
% load([file_name '\foreground_iter1_jointThresholding2\score_map_' num2str(tt+1)]);
% [score_max, score_max_ind] = max(score_map,[],4);
% score_max_ind(score_max==0) = 0;

% get relationship
% foreground = imregionalmax(score_max,conn);
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
tifwrite(uint8(result*255),[file_name '\' method_name '\' num2str(smoothness) '_' num2str(z_thre) '\' num2str(tt+1)]);
end