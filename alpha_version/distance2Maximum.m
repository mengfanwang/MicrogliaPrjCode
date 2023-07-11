clc;clear;close all;
dbstop if error
% distance_map_max distance map after nms
% tip1 distance_map_max with smoothing
% tip2 local maximum region of tip

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\redetection_info_v2\tipDetection';
% mkdir([file_name 'tip_distance\distance_map_max']);
% mkdir([file_name 'tip_distance\ball']);
mkdir(fullfile(file_name, 'distance_map_max'));
mkdir(fullfile(file_name, 'fore_ball'));
mkdir(fullfile(file_name, 'data_ball'));
t = 33;

xCoord = cell(t,1);
yCoord = cell(t,1);
zCoord = cell(t,1);
tip_score = cell(t,1);
com_label = cell(t,1);
for tt = 1:t
    tt
ind = num2str(tt+1);
% ind = ind(2:4);
% load([file_name 'tip_distance\' ind '.mat']);
load(fullfile(file_name, [ind '.mat']));
[x_size, y_size, z_size, step] = size(distance_map);
step = 11:20;
distance_map_max = zeros(x_size, y_size, z_size, length(step));
for ii = 1:length(step)
    distance_map_max(:,:,:,ii) = distance_map(:,:,:,step(ii))/step(ii);
end
distance_map_max = max(distance_map_max,[],4);
% writeColormap(distance_map_max, [file_name 'tip_distance\distance_map_max\'  ind]);
writeColormap(distance_map_max, fullfile(file_name, 'distance_map_max', ind));

% load data
distance_map = distance_map_max;
fore_all = ~isnan(distance_map);
com = bwconncomp(fore_all);
com_size = cellfun(@length,com.PixelIdxList);
distance_max = nan(size(fore_all));
distance_max_ind = [];

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

tic;
tip1 = nan(size(fore_all));
for ii = 1:com.NumObjects
    
    distance_temp = zeros(size(fore_all));
    distance_temp2 = zeros(size(fore_all));
    distance_temp(com.PixelIdxList{ii}) = distance_map(com.PixelIdxList{ii});
    [x_ind, y_ind ,z_ind] = ind2sub(size(fore_all),com.PixelIdxList{ii});
    % use mean filter to denoise
    for jj = 1:length(com.PixelIdxList{ii})
        if distance_map(com.PixelIdxList{ii}(jj)) > 0
        distance_win = distance_temp(max(x_ind(jj) - 2, 1):min(x_ind(jj) + 2, x_size),...
            max(y_ind(jj) - 2, 1):min(y_ind(jj) + 2, y_size), max(z_ind(jj) - 2, 1):min(z_ind(jj) + 2, z_size));
        pix_geodesic = getGeodesicPixel(distance_win, [x_ind(jj) y_ind(jj) z_ind(jj)], 2);
        distance_temp2(x_ind(jj), y_ind(jj), z_ind(jj)) = mean(distance_win(pix_geodesic));
        tip1(x_ind(jj), y_ind(jj), z_ind(jj)) = mean(distance_win(pix_geodesic));
        else
        distance_temp2(x_ind(jj), y_ind(jj), z_ind(jj)) = 0;
        tip1(x_ind(jj), y_ind(jj), z_ind(jj)) = 0;
        end
    end
    distance_temp = distance_temp2;
    for jj = 1:length(com.PixelIdxList{ii})
        if distance_map(com.PixelIdxList{ii}(jj)) <= 0
            continue
        end
        distance_win = distance_temp(max(x_ind(jj) - 5, 1):min(x_ind(jj) + 5, x_size),...
            max(y_ind(jj) - 5, 1):min(y_ind(jj) + 5, y_size), max(z_ind(jj) - 5, 1):min(z_ind(jj) + 5, z_size));
        pix_geodesic = getGeodesicPixel(distance_win, [x_ind(jj) y_ind(jj) z_ind(jj)], 5);
        if max(distance_win(pix_geodesic)) == distance_temp(x_ind(jj), y_ind(jj), z_ind(jj)) && max(distance_win(pix_geodesic)) > 0.5
            distance_max(x_ind(jj), y_ind(jj), z_ind(jj)) = distance_temp(x_ind(jj), y_ind(jj), z_ind(jj));
            % region growing under certain std
            max_std = std(distance_win(distance_win>0));
            thr = distance_temp(x_ind(jj), y_ind(jj), z_ind(jj)) - max_std;
            candidate_list = [x_ind(jj) y_ind(jj) z_ind(jj)];
            while ~isempty(candidate_list)
                candidate = candidate_list(1,:);
                candidate_list(1,:) = [];
                for dd = 1:26
                    x_neighbor = min(max(candidate(1) + x_direction(dd),1),x_size); 
                    y_neighbor = min(max(candidate(2) + y_direction(dd),1),y_size); 
                    z_neighbor = min(max(candidate(3) + z_direction(dd),1),z_size); 
                    if isnan(distance_max(x_neighbor, y_neighbor, z_neighbor)) && ...
                            distance_temp(x_neighbor, y_neighbor, z_neighbor) > thr
                        distance_max(x_neighbor, y_neighbor, z_neighbor) = distance_temp(x_neighbor, y_neighbor, z_neighbor);
                        candidate_list = [candidate_list; x_neighbor y_neighbor z_neighbor];
                    end
                end
            end
        end
    end
end
toc


% write image
% writeColormap(tip1, [file_name '\tip_distance\']);
% x = distance_max(distance_max>0);
% fore_tip = zeros(size(fore_all));
% fore_tip(distance_max>0) = x;
% fore_tip = imdilate(fore_tip,ones(3,3,3));
% im = zeros(x_size, y_size, 3, z_size);
% for zz = 1:z_size
%     im(:,:,1,zz) = fore_all(:,:,zz)*128;
%     im(:,:,2,zz) = fore_all(:,:,zz)*128;
%     im(:,:,3,zz) = fore_all(:,:,zz)*128;
% end
% fore_tip = ceil(fore_tip*63/max(fore_tip(:))-1e-8) + 1;
% color = colormap;
% for ii = 1:x_size
%     for jj = 1:y_size
%         for kk = 1:z_size
%             if fore_tip(ii,jj,kk) > 1
%                 im(ii,jj,1,kk) = color(fore_tip(ii,jj,kk),1)*255;
%                 im(ii,jj,2,kk) = color(fore_tip(ii,jj,kk),2)*255;
%                 im(ii,jj,3,kk) = color(fore_tip(ii,jj,kk),3)*255;
%             end
%         end
%     end
% end
% tifwrite(uint8(im), [file_name '\tip_distance\tip2']);


% x = distance_max(distance_max>0);
fore_tip = zeros(size(fore_all));
% fore_tip(distance_max>0) = x - 0.5;

% foreground = load([file_name '\foreground_c2\' ind '.mat']);
foreground = tifread(fullfile('D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\redetection_info_v2\foreground', [ind '.tif']));
foreground = logical(foreground);
conn_label = bwlabeln(foreground);
com = bwconncomp(distance_max>0);
for ii = 1:com.NumObjects
    [x_ind, y_ind ,z_ind] = ind2sub(size(fore_all),com.PixelIdxList{ii});
    distance_temp = distance_max(com.PixelIdxList{ii});
    distance_norm = distance_temp/sum(distance_temp);
    x_ind = round(sum(x_ind.*distance_norm));
    y_ind = round(sum(y_ind.*distance_norm));
    z_ind = round(sum(z_ind.*distance_norm));
    fore_tip(x_ind, y_ind, z_ind) = mean(distance_temp);
    label_ind = unique(conn_label(com.PixelIdxList{ii}));
    if length(label_ind) ~= 1
        error('Some tip are in different components.');
    end
    xCoord{tt} = [xCoord{tt}; x_ind];
    yCoord{tt} = [yCoord{tt}; y_ind];
    zCoord{tt} = [zCoord{tt}; z_ind];
    tip_score{tt} = [tip_score{tt}; mean(distance_temp)];
    com_label{tt} = [com_label{tt}; label_ind];
end
unique_flag = zeros(length(com_label{tt}),1);
for ii = 1:length(com_label{tt})
    if sum(com_label{tt} == com_label{tt}(ii)) == 1
        unique_flag(ii) = 1;
    end
end
com_label{tt} = [com_label{tt} unique_flag];

fore_tip = imdilate(fore_tip,strel('sphere',3));
im_fore = zeros(x_size, y_size, 3, z_size);
ind2 = num2str(1001+tt);
ind2 = ind2(2:4);
data = tifread(fullfile('D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\time_domain_data_after_VST', [ind2 '.tif']));
data = double(data);
im_data = zeros(x_size, y_size, 3, z_size);
for zz = 1:z_size
    im_fore(:,:,1,zz) = fore_all(:,:,zz)*128;
    im_fore(:,:,2,zz) = fore_all(:,:,zz)*128;
    im_fore(:,:,3,zz) = fore_all(:,:,zz)*128;
    im_data(:,:,1,zz) = data(:,:,zz);
    im_data(:,:,2,zz) = data(:,:,zz);
    im_data(:,:,3,zz) = data(:,:,zz);
end
fore_tip = ceil(fore_tip*63/max(fore_tip(:))-1e-8) + 1;
color = colormap;
for ii = 1:x_size
    for jj = 1:y_size
        for kk = 1:z_size
            if fore_tip(ii,jj,kk) > 1
                im_fore(ii,jj,1,kk) = color(fore_tip(ii,jj,kk),1)*255;
                im_fore(ii,jj,2,kk) = color(fore_tip(ii,jj,kk),2)*255;
                im_fore(ii,jj,3,kk) = color(fore_tip(ii,jj,kk),3)*255;
                im_data(ii,jj,1,kk) = color(fore_tip(ii,jj,kk),1)*255;
                im_data(ii,jj,2,kk) = color(fore_tip(ii,jj,kk),2)*255;
                im_data(ii,jj,3,kk) = color(fore_tip(ii,jj,kk),3)*255;
            end
        end
    end
end
% tifwrite(uint8(im), [file_name 'tip_distance\ball\'  ind]);
tifwrite(uint8(im_fore), fullfile(file_name, 'fore_ball', ind));
tifwrite(uint8(im_data), fullfile(file_name, 'data_ball', ind));
end
% save([file_name 'tip_info.mat'], 'xCoord', 'yCoord', 'zCoord', 'tip_score', 'com_label');
save(fullfile(file_name, 'tip_info.mat'), 'xCoord', 'yCoord', 'zCoord', 'tip_score', 'com_label');

function pix_geodesic = getGeodesicPixel(distance_win, idx, step)
    center_idx = min(idx,step+1);
    center_idx = sub2ind(size(distance_win), center_idx(1), center_idx(2), center_idx(3));
    distance_geodesic = bwdistgeodesic(distance_win>0,center_idx,'quasi-euclidean');
    pix_geodesic = distance_geodesic < step;
end

function writeColormap(distance_map, filename)
distance_map(distance_map>0.8) = 0.8;
distance_map = ceil(distance_map*63/max(distance_map(:))-1e-8) + 1;
color = colormap;
[x,y,z] = size(distance_map);
map = zeros(x,y,3,z);
for ii = 1:x
    for jj = 1:y
        for kk = 1:z
            if ~isnan(distance_map(ii,jj,kk))
                map(ii,jj,1,kk) = color(distance_map(ii,jj,kk),1);
                map(ii,jj,2,kk) = color(distance_map(ii,jj,kk),2);
                map(ii,jj,3,kk) = color(distance_map(ii,jj,kk),3);
            end
        end
    end
end
tifwrite(uint8(map*255),filename);
end