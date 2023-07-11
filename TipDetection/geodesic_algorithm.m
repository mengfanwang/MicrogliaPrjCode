clc;clear;close all;
x_size = 512; y_size = 511; z_size = 61;
% the last column is deleted becasue of the data property
% in the original data, the last column is almost empty
dbstop if error 
% addpath C:\Users\mengf\Desktop\geodesic\geodesic

% for_all bound_all node_num bound_num win_size

win_size = 15;
neighbor = cat(3, [0 0 0; 0 1 0; 0 0 0],...
                  [0 1 0; 1 1 1; 0 1 0],...
                  [0 0 0; 0 1 0; 0 0 0]);
% calculate distance
step = (win_size-1)/2;
center = [(win_size+1)/2 (win_size+1)/2 (win_size+1)/2];
center_idx = sub2ind([win_size win_size win_size],center(1),center(2),center(3));
              
load fore_all.mat
fore_all = double(fore_all);
bound_temp = fore_all - imerode(fore_all, neighbor);
bound_all = zeros(size(fore_all));
bound_all(step+1:end-step,step+1:end-step,step+1:end-step) = bound_temp(step+1:end-step,step+1:end-step,step+1:end-step);

fore_ind = find(fore_all);
node_num = length(fore_ind);
fore_all(fore_ind) = 1:node_num;
bound_ind = find(bound_all);
bound_num = length(bound_ind);
% dlmwrite('C:\Users\mengf\Desktop\geodesic\geodesic_v2\fore_all.dat',node_num,'precision','%d');
% for ii = 1:z_size
%     dlmwrite('C:\Users\mengf\Desktop\geodesic\geodesic_v2\fore_all.dat',fore_all(:,:,ii),'-append','delimiter',' ','precision','%d');
% end
% dlmwrite('C:\Users\mengf\Desktop\geodesic\geodesic_v2\bound_all.dat',bound_num,'precision','%d');
% for ii = 1:z_size
%     dlmwrite('C:\Users\mengf\Desktop\geodesic\geodesic_v2\bound_all.dat',bound_all(:,:,ii),'-append','delimiter',' ','precision','%d');
% end
[x_ind,y_ind,z_ind] = ind2sub(size(bound_all),find(bound_all));
bound_list = [x_ind y_ind z_ind fore_all(bound_ind)];
% node(1:node_num) = struct('edge',zeros(26,2),'edge_num',0,'dist',inf,'flag',0);
node_edge = zeros(node_num,26);
node_capacity = zeros(node_num,26);
node_edge_num = zeros(node_num,1);
node_dist = inf(node_num,1);
node_flag = zeros(node_num,1);

direction = zeros(26,4);
direction_ind = 0;
for xx = -1:1
    for yy = -1:1
        for zz = -1:1
            if abs(xx) + abs(yy) + abs(zz) > 0
                direction_ind = direction_ind + 1;
                direction(direction_ind,:) = [xx yy zz abs(xx) + abs(yy) + abs(zz)];
            end
        end
    end
end
c = [1 sqrt(2) sqrt(3)];
tic;
for x = 1:x_size
    x
    for y = 1:y_size
        for z = 1:z_size
            if fore_all(x,y,z) > 0
                node_s = fore_all(x,y,z);
                for dd = 1:26
                    xx = x+direction(dd,1);
                    yy = y+direction(dd,2);
                    zz = z+direction(dd,3);
                    if xx>=1 && xx<=x_size && yy>=1 && yy<=y_size && zz>=1 && zz<=z_size
                        if fore_all(xx,yy,zz) > 0
                            node_t = fore_all(xx,yy,zz);
%                             node(node_s).edge_num = node(node_s).edge_num + 1;
%                             node(node_s).edge(node(node_s).edge_num,:) = [node_t direction(dd,4)];
                            node_edge_num(node_s) = node_edge_num(node_s) + 1;
                            node_edge(node_s,node_edge_num(node_s)) = node_t;
                            node_capacity(node_s,node_edge_num(node_s)) = direction(dd,4);
                        end
                    end
                end
            end
        end
    end
end
% toc

% load node_struct.mat
% load edge_list.mat
% edge_list4 = zeros(10000,1);
% for ii = 1:10000
%     edge_list4(ii) = fore_all(edge_list(ii,1),edge_list(ii,2),edge_list(ii,3));
% end
% bound_list = [edge_list edge_list4];

tic;
time = 0;
% for jj = 1:10000
jj = 31342;
win = fore_all(bound_list(jj,1)-step:bound_list(jj,1)+step,bound_list(jj,2)-step:bound_list(jj,2)+step,bound_list(jj,3)-step:bound_list(jj,3)+step);
dist_map2 = bwdistgeodesic(logical(win),center_idx,'quasi-euclidean');

node_ind = find(win);
node_list = win(node_ind);
for ii = 1:length(node_ind)
%     node(node_list(ii)).flag = 1;
    node_flag(node_list(ii)) = 1;
end

current_time = 0;
edge_queue = zeros(26*numel(win),4);
node_queue = zeros(numel(win),1);
edgeind = 0;

node_s = win(center(1),center(2),center(3));
node_dist(node_s) = 0;
for ee = 1:node_edge_num(node_s)
    edgeind = edgeind + 1;
    edge_queue(edgeind,1:2) = [node_s node_edge(node_s,ee)];
    edge_queue(edgeind,4) = c(node_capacity(node_s,ee));
end
while edgeind > 0
current_time = current_time + 1;
edgeind_new = 0;
nodeind = 0;
for ii = 1:edgeind
    edge_queue(ii,3) = edge_queue(ii,3) + 1;
    if edge_queue(ii,3) < edge_queue(ii,4)
        edgeind_new = edgeind_new + 1;
        edge_queue(edgeind_new,:) = edge_queue(ii,:);
    else
        fire_time = current_time - (edge_queue(ii,3) - edge_queue(ii,4));
        node_t = edge_queue(ii,2);
        if node_dist(node_t) == inf
            nodeind = nodeind + 1;
            node_queue(nodeind) = node_t;
            node_dist(node_t) = fire_time;
        elseif node_dist(node_t) > fire_time
            node_dist(node_t) = fire_time;
        end
        disp([num2str(node_t) ' ' num2str(fire_time)]);
    end
end
% tic;
for ii = 1:nodeind
    node_s = node_queue(ii);
    for ee = 1:node_edge_num(node_s)
        node_t = node_edge(node_s,ee);
%         a = a + 1;
        if node_flag(node_t) == 1 && node_dist(node_t) == inf        
            edgeind_new = edgeind_new + 1;
            edge_queue(edgeind_new,:) =[node_s node_t current_time-node_dist(node_s) c(node_capacity(node_s,ee))];
        end
    end
end
% time = time + toc;
edgeind = edgeind_new;
end

dist_map = nan(size(win));
for ii = 1:length(node_ind)
    dist_map(node_ind(ii)) = node_dist(node_list(ii));
    node_flag(node_list(ii)) = 0;
    node_dist(node_list(ii)) = inf;
end
 

% if sum(dist_map(:)<100) - sum(dist_map2(:)<100) ~= 0
%     error('da');
% end
% if max(dist_map(:) - dist_map2(:))>1e-4
%     error('daa');
% end
% end
toc