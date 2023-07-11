clc;clear;close all;
dbstop if error

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';
mkdir([file_name '\tip_connection']);
load([file_name 'data_c2.mat']);
data_all = data2;
load([file_name 'tip_info.mat']);
[x, y, z, t] = size(data_all);

%%% parameter setting %%%
radius = 15;
exter = 5; % consider bending branches
max_path = 0.2;


for tt = 1:t
    tt
data = data_all(:,:,:,tt);
tip_num = length(xCoord{tt});
mindist_matrix = inf(size(data));
connection_num = 0;
d_ = nan(tip_num^2,1);
for ii = 1:tip_num
    for jj = ii+1:tip_num
        tip1 = [xCoord{tt}(ii) yCoord{tt}(ii) zCoord{tt}(ii)];
        tip2 = [xCoord{tt}(jj) yCoord{tt}(jj) zCoord{tt}(jj)];
        if norm(tip1 - tip2) < radius
            tip1 = round(tip1);
            tip2 = round(tip2);
            if com_label{tt}(ii,1) ~= com_label{tt}(jj,1)
                x_min = max(min(tip1(1), tip2(1)) - exter, 1);
                y_min = max(min(tip1(2), tip2(2)) - exter, 1);
                z_min = max(min(tip1(3), tip2(3)) - exter, 1);
                x_max = min(max(tip1(1), tip2(1)) + exter, x);
                y_max = min(max(tip1(2), tip2(2)) + exter, y);
                z_max = min(max(tip1(3), tip2(3)) + exter, z);
                limits = [x_min y_min z_min x_max y_max z_max];
                [short_path, d] = findPath(tip1, tip2, data, limits);
                d_(jj + tip_num*ii) = d;
                if d < max_path
                    connection_num = connection_num + 1;
                    connection_info(connection_num).source = [xCoord{tt}(ii) yCoord{tt}(ii) zCoord{tt}(ii)];
                    connection_info(connection_num).sink   = [xCoord{tt}(jj) yCoord{tt}(jj) zCoord{tt}(jj)];
                    connection_info(connection_num).distance = d;
                    connection_info(connection_num).shortestpath = short_path;    
                    connection_info(connection_num).source_ind = ii;
                    connection_info(connection_num).sink_ind = jj;
                    
                    mindist_matrix(tip1(1),tip1(2),tip1(3)) = min( mindist_matrix(tip1(1),tip1(2),tip1(3)), d);
                    mindist_matrix(tip2(1),tip2(2),tip2(3)) = min( mindist_matrix(tip2(1),tip2(2),tip2(3)), d);
                end
            end
        end
    end
end
% remove duplicate connections
remove_list = [];
for ii = 1:length(connection_info)
    tip1 = round(connection_info(ii).source);
    tip2 = round(connection_info(ii).sink);
    if mindist_matrix(tip1(1),tip1(2),tip1(3)) ~= connection_info(ii).distance || mindist_matrix(tip2(1),tip2(2),tip2(3)) ~= connection_info(ii).distance
        remove_list = [remove_list ii];
    end
end
connection_info(remove_list) = [];
% keep useful tip (not connected/connected but single)
remove_list = zeros(length(xCoord{tt}),1);
for ii = 1:length(connection_info)
    if com_label{tt}(connection_info(ii).source_ind,2) == 0
        remove_list(connection_info(ii).source_ind) = 1;
    end
    if com_label{tt}(connection_info(ii).sink_ind,2) == 0
        remove_list(connection_info(ii).sink_ind) = 1;
    end
end
remove_list = logical(remove_list);
xCoord{tt}(remove_list) = [];
yCoord{tt}(remove_list) = [];
zCoord{tt}(remove_list) = [];
tip_score{tt}(remove_list) = [];
com_label{tt}(remove_list,:) = [];

% plot part
ind = num2str(1000+tt);
ind = ind(2:4);
foreground = load([file_name 'foreground_c2\' ind '.mat']);
foreground = foreground.foreground;
% foreground = data_all(:,:,:,tt)/128;
% foreground = tifread([file_name '\stabilization\' ind '.tif'])/128;

% plot connection (line)
im = zeros(x, y, 3, z);
for ii = 1:length(connection_info)
%     s = round(connection_info(ii).source);
%     t = round(connection_info(ii).sink);
%     [im, label_pix] = draw3Dline(im, s, t, 1.5, [1 1 1]*255);
    shortest_path = connection_info(ii).shortestpath;
    for jj = 1:size(shortest_path,1)
        im(shortest_path(jj,1),shortest_path(jj,2),1,shortest_path(jj,3)) = 255;
        im(shortest_path(jj,1),shortest_path(jj,2),2,shortest_path(jj,3)) = 255;
        im(shortest_path(jj,1),shortest_path(jj,2),3,shortest_path(jj,3)) = 255;
    end
end
conn = cat(3, [0 1 0; 1 1 1; 0 1 0],...
              [1 1 1; 1 1 1; 1 1 1],...
              [0 1 0; 1 1 1; 0 1 0]);
im = imdilate(im,conn);
for zz = 1:z
    im(:,:,1,zz) = im(:,:,1,zz) + foreground(:,:,zz)*128;
    im(:,:,2,zz) = im(:,:,2,zz) + foreground(:,:,zz)*128;
    im(:,:,3,zz) = im(:,:,3,zz) + foreground(:,:,zz)*128;
end
% tifwrite(uint8(im), [file_name '\tip_connection\'  ind]);


fore_tip = zeros(size(foreground));
for ii = 1:length(xCoord{tt})
    x_ind = xCoord{tt}(ii);
    y_ind = yCoord{tt}(ii);
    z_ind = zCoord{tt}(ii);
    fore_tip(x_ind, y_ind, z_ind) = tip_score{tt}(ii);
end
fore_tip = imdilate(fore_tip,strel('sphere',3));
% im = zeros(x, y, 3, z);
% for zz = 1:z
%     im(:,:,1,zz) = foreground(:,:,zz)*128;
%     im(:,:,2,zz) = foreground(:,:,zz)*128;
%     im(:,:,3,zz) = foreground(:,:,zz)*128;
% end
fore_tip = ceil(fore_tip*63/max(fore_tip(:))-1e-8) + 1;
color = colormap;
for ii = 1:x
    for jj = 1:y
        for kk = 1:z
            if fore_tip(ii,jj,kk) > 1
                im(ii,jj,1,kk) = color(fore_tip(ii,jj,kk),1)*255;
                im(ii,jj,2,kk) = color(fore_tip(ii,jj,kk),2)*255;
                im(ii,jj,3,kk) = color(fore_tip(ii,jj,kk),3)*255;
            end
        end
    end
end
tifwrite(uint8(im), [file_name '\tip_connection\'  ind]);
end
save([file_name 'tip_info_reconnect.mat'], 'xCoord', 'yCoord', 'zCoord', 'tip_score', 'com_label');

function [short_path, d] = findPath(tip1, tip2, data, limits)
    x_direction = zeros(1,26);
    y_direction = zeros(1,26);
    z_direction = zeros(1,26);
    ind_dir = 0;
    for xx = -1:1
        for yy = -1:1
            for zz = -1:1
                if abs(xx) + abs(yy) + abs(zz) > 0
                    ind_dir = ind_dir + 1;
                    x_direction(ind_dir) = xx;
                    y_direction(ind_dir) = yy;
                    z_direction(ind_dir) = zz;
                end
            end
        end
    end
    [x, y, z] = size(data);
    x_min = limits(1); y_min = limits(2); z_min = limits(3);
    x_max = limits(4); y_max = limits(5); z_max = limits(6);
    x_length = x_max - x_min + 1;
    y_length = y_max - y_min + 1;
    z_length = z_max - z_min + 1;
    node_table = zeros(x_length*y_length*z_length,3);
    node_matrix = zeros(size(data));
    for xx = x_min:x_max
        for yy = y_min:y_max
            for zz = z_min:z_max
                node_ind = (xx - x_min) + (yy - y_min)*x_length + (zz - z_min)*x_length*y_length + 1;
                node_table(node_ind,:) = [xx yy zz];
                node_matrix(xx,yy,zz) = node_ind;
            end
        end
    end
    % build graph
    node_num = node_ind;
    edge_num = (x_length-2)*(y_length-2)*(z_length-2)*26 + 2*( (x_length-2)*(y_length-2) + (x_length-2)*(z_length-2) + ...
                (y_length-2)*(z_length-2) )*17 + 4*( (x_length-2)+(y_length-2)+(z_length-2) )*11 + 8*7;
    s = zeros(edge_num,1);
    t = zeros(edge_num,1);
    w = zeros(edge_num,1);
    edge_ind = 0;
    for ii = 1:node_num
        node_ind = node_table(ii,:);
        for dd = 1:26
            x_neighbor = node_ind(1) + x_direction(dd);
            y_neighbor = node_ind(2) + y_direction(dd);
            z_neighbor = node_ind(3) + z_direction(dd);
            if x_neighbor >= 1 && x_neighbor <= x && y_neighbor >= 1 && y_neighbor <= y && z_neighbor >= 1 && z_neighbor <= z
                if node_matrix(x_neighbor, y_neighbor, z_neighbor) > 0
                    edge_ind = edge_ind + 1;
                    s(edge_ind) = ii;
                    t(edge_ind) = node_matrix(x_neighbor, y_neighbor, z_neighbor);
                    w(edge_ind) = 1/(data(x_neighbor, y_neighbor, z_neighbor)+1);
                end
            end
        end
    end
    g = digraph(s,t,w);
    [p,d] = shortestpath(g, node_matrix(tip1(1), tip1(2), tip1(3)), node_matrix(tip2(1), tip2(2), tip2(3)));
    short_path = zeros(length(p), 3);
    for ii = 1:length(p)
        short_path(ii,:) = [node_table(p(ii),1) node_table(p(ii),2) node_table(p(ii),3)];
    end
end