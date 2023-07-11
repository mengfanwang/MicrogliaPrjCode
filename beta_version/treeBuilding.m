function [result, tip_info, bran_info, tip2bran, Node, G] = treeBuilding(foreground,soma)
[x_size, y_size, z_size] = size(foreground);
dis = bwdistgeodesic(foreground,soma,'quasi-euclidean');
dis(isnan(dis)) = 0;
dis(dis==Inf) = 0;
foreground = logical(dis) + soma;

dis_floor = ceil(dis/2);
dis_list = unique(dis_floor);
dis_list = dis_list(2:end); % remove 0

for ii = 1:length(dis_list)
    dis_layer(ii) = bwconncomp(dis_floor == dis_list(ii),26);
end

node_map = double(soma);
node_label = soma;
node_ind = 1;
Node(node_ind).node = node_ind;
Node(node_ind).parent = [];
Node(node_ind).child = [];
Node(node_ind).voxel = find(node_map == node_ind);
Node(node_ind).distance = 0;
[x_ind,y_ind,z_ind] = ind2sub(size(foreground),Node(node_ind).voxel);
Node(node_ind).x = mean(x_ind);
Node(node_ind).y = mean(y_ind);
Node(node_ind).z = mean(z_ind);
Node(node_ind).color = 1;
G = digraph;

% tic;
candidate_list = find((imdilate(node_label,ones(3,3,3)) - node_label)==1);
candidate_list(dis_floor(candidate_list)==0) = [];
while ~isempty(candidate_list)
    while ~isempty(candidate_list)
        candidate = candidate_list(1);
        candidate_list(1) = [];
        % if detected we don't take option
        if node_map(candidate) ~= 0
            continue
        end
        candidate_dis = dis_floor(candidate);
        % build a new node
        for ii = 1:dis_layer(candidate_dis).NumObjects
            if ismember(candidate,dis_layer(candidate_dis).PixelIdxList{ii})
                node_ind = node_ind + 1;
                Node(node_ind).node = node_ind;
                Node(node_ind).parent = [];
                Node(node_ind).child = [];
                Node(node_ind).voxel = dis_layer(candidate_dis).PixelIdxList{ii};
                Node(node_ind).distance = candidate_dis;
                [x_ind,y_ind,z_ind] = ind2sub(size(foreground),Node(node_ind).voxel);
                Node(node_ind).x = mean(x_ind);
                Node(node_ind).y = mean(y_ind);
                Node(node_ind).z = mean(z_ind);
                node_map(Node(node_ind).voxel) = node_ind;
                node_label(Node(node_ind).voxel) = 1;
%                 % delete the component
%                 dis_layer(candidate_dis).PixelIdxList(ii) = [];
%                 dis_layer(candidate_dis).NumObjects = dis_layer(candidate_dis).NumObjects - 1;
            end
        end
        % find the parent
        candidate_region = zeros(size(foreground));
        candidate_region(Node(node_ind).voxel) = 1;
        parent_cand_list = find((imdilate(candidate_region,ones(3,3,3)) - candidate_region)==1);
        for ii = 1:length(parent_cand_list)
            parent_cand = parent_cand_list(ii);
            if node_map(parent_cand) ~= 0 && dis_floor(parent_cand) < candidate_dis &&...
                ~ismember(node_ind,Node(node_map(parent_cand)).child)
                Node(node_map(parent_cand)).child = [Node(node_map(parent_cand)).child node_ind];
                Node(node_ind).parent = [Node(node_ind).parent node_map(parent_cand)];
                Node(node_ind).color = 1 - Node(node_map(parent_cand)).color;
                G = addedge(G,node_map(parent_cand),node_ind); 
            end
        end
    end
    candidate_list = find((imdilate(node_label,ones(3,3,3)) - node_label)==1);
    candidate_list(dis_floor(candidate_list)==0) = [];
end
% toc;

% find leaves
leaves = [];
nLabel = zeros(1,length(Node));
x_ind = zeros(1,length(Node));
y_ind = zeros(1,length(Node));
z_ind = zeros(1,length(Node));
for ii = 1:length(Node)
    nLabel(ii) = Node(ii).distance;
    x_ind(ii) = Node(ii).x;
    y_ind(ii) = Node(ii).y;
    z_ind(ii) = Node(ii).z;
    if isempty(Node(ii).child)
        leaves = [leaves ii];
        Node(ii).branch_len = 0;
    end
end

% build color
result_red = zeros(size(foreground));
result_blue = zeros(size(foreground));
for ii = 1:length(Node) 
    if Node(ii).color == 1
        result_red(Node(ii).voxel) = 1;
    elseif Node(ii).color == 0
        result_blue(Node(ii).voxel) = 1;
    end
end

% find branches from leaves
for ii = 1:length(Node)
    if isempty(Node(ii).branch_len)
        Node(ii).branch_len = zeros(size(Node(ii).child));
    end
end
for ii = length(Node):-1:1
    for jj = 1:length(Node(ii).parent)
        Node(Node(ii).parent(jj)).branch_len(Node(Node(ii).parent(jj)).child == ii) = max(Node(ii).branch_len) + 1;
    end
end
                
        


% figure(1);
% h = plot(G,'XData',x_ind,'YData',y_ind,'ZData',z_ind);%,'NodeLabel',nLabel);%,'Layout','force');
% highlight(h,leaves,'NodeColor','r');
% 
nLabel2 = zeros(1,length(Node));
for ii = 1:length(Node)
    nLabel2(ii) = max(Node(ii).branch_len);
end
% 
% figure(2);
% h = plot(G,'NodeLabel',nLabel2,'Layout','force');
% highlight(h,leaves,'NodeColor','r');   

% remove small branch
tip_label = ones(1,length(leaves));
for ii = 1:length(leaves)
%     % find parent
%     child_ind = leaves(ii);
%     branch_ind = Node(child_ind).parent;
%     if length(branch_ind) > 1
%         continue
%     end
%     while 1
%         if length(Node(branch_ind).child) > 1
%             break
%         end
%         if length(Node(branch_ind).parent) == 1
%             child_ind = branch_ind;
%             branch_ind = Node(branch_ind).parent;
%         else
%             break
%         end
%     end
%     if (max(Node(branch_ind).branch_len) > 3 &&...
%             Node(branch_ind).branch_len(Node(branch_ind).child == child_ind) < 3) || ...
%        (max(Node(branch_ind).branch_len) > 2 &&...    
%             Node(branch_ind).branch_len(Node(branch_ind).child == child_ind) < 2)
%         tip_label(ii) = 0;
%     end
    
    % remove small distance
    if Node(leaves(ii)).distance < 10
        tip_label(ii) = 0;
    end
end
leaves(tip_label==0) = [];

% tip to branch length
tip2bran = zeros(1,length(leaves));
for ii = 1:length(leaves)
    node_ind = leaves(ii);
    while 1
        tip2bran(ii) = tip2bran(ii)+1;
        if length(Node(node_ind).parent) > 1 || length(Node(Node(node_ind).parent).child) > 1
            break
        end
        node_ind = Node(node_ind).parent;
    end
end
        

% add branch
branch = [];
for ii = 1:length(Node)
    if sum(Node(ii).branch_len) >= 10 && length(Node(ii).branch_len) > 2 && Node(ii).distance > 10
        branch = [branch ii];
    end
end

% nLabel3 = zeros(1,length(Node));
% for ii = 1:length(Node)
%     nLabel3(ii) = Node(ii).color;
% end
% nLabel3 = find(nLabel3);
% 
% figure(3); 
% h1 = plot(G,'XData',x_ind,'YData',y_ind,'ZData',z_ind);%,'NodeLabel',nLabel);%,'Layout','force');
% highlight(h1,nLabel3,'NodeColor','r');
% highlight(h1,leaves,'NodeColor','y');
% highlight(h1,branch,'NodeColor','g');
% set(gca,'Color','k');
% 
% figure(4);
% h2 = plot(G,'NodeLabel',nLabel2,'Layout','force');
% highlight(h2,leaves,'NodeColor','r');
  
% tip & branch
tip = zeros(size(foreground));
tip_info = zeros(length(leaves),3);
for ii = 1:length(leaves)
%     tip(Node(leaves(ii)).voxel) = 1;
    tip(round(Node(leaves(ii)).x),round(Node(leaves(ii)).y),round(Node(leaves(ii)).z)) = 1;
    tip_info(ii,:) = [Node(leaves(ii)).x Node(leaves(ii)).y Node(leaves(ii)).z];
end
tip = imdilate(tip,strel('sphere',2));
bran = zeros(size(foreground));
bran_info = zeros(length(branch),3);
for ii = 1:length(branch)
    bran(round(Node(branch(ii)).x),round(Node(branch(ii)).y),round(Node(branch(ii)).z)) = 1;
    bran_info(ii,:) = [Node(branch(ii)).x Node(branch(ii)).y Node(branch(ii)).z];
end
bran = imdilate(bran,strel('sphere',3));

result = zeros(x_size,y_size,3,z_size);
for zz = 1:z_size
    result(:,:,1,zz) = result_red(:,:,zz)*255;
    result(:,:,3,zz) = result_blue(:,:,zz)*255;
%     result(:,:,1,zz) = max(tip(:,:,zz)*255,result_red(:,:,zz)*255);
%     result(:,:,2,zz) = tip(:,:,zz)*255;
%     result(:,:,3,zz) = result_blue(:,:,zz)*255;
end

% [x_ind2,y_ind2,z_ind2] = ind2sub(size(foreground),find(foreground));
% x_max = max(x_ind2); x_min = min(x_ind2);
% y_max = max(y_ind2); y_min = min(y_ind2);
% z_max = max(z_ind2); z_min = min(z_ind2);
% tifwrite(uint8(result(max(x_min-5,1):min(x_max+5,x_size),max(y_min-5,1):min(y_max+5,y_size),:,...
%     max(z_min-5,1):min(z_max+5,z_size))),'result_remove');