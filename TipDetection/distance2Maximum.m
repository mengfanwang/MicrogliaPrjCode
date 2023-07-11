clc;clear;close all;
% keep the local maximum distances only

%kernel_size = 3;
x_size = 512; y_size = 511; z_size = 61;
% load data
load distance_map_max.mat
distance_map = distance_map_max;
load fore_all.mat
com = bwconncomp(fore_all);
com_size = cellfun(@length,com.PixelIdxList);
distance_max = nan(size(fore_all));
distance_max_ind = [];
fore_temp = zeros(size(fore_all));

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

for ii = 1:com.NumObjects
    if com_size(ii) > 2000 
        ii
        distance_temp = zeros(size(fore_all));
        distance_temp2 = zeros(size(fore_all));
        distance_temp(com.PixelIdxList{ii}) = distance_map(com.PixelIdxList{ii});
        fore_temp(com.PixelIdxList{ii}) = fore_all(com.PixelIdxList{ii});
        [x_ind, y_ind ,z_ind] = ind2sub(size(fore_all),com.PixelIdxList{ii});
        % use mean filter to denoise
        for jj = 1:length(com.PixelIdxList{ii})
            distance_win = distance_temp(max(x_ind(jj) - 2, 1):min(x_ind(jj) + 2, x_size),...
                max(y_ind(jj) - 2, 1):min(y_ind(jj) + 2, y_size), max(z_ind(jj) - 2, 1):min(z_ind(jj) + 2, z_size));
            distance_temp2(x_ind(jj), y_ind(jj), z_ind(jj)) = mean(distance_win(distance_win>0));
        end
        distance_temp = distance_temp2;
        for jj = 1:length(com.PixelIdxList{ii})
            distance_win = distance_temp(max(x_ind(jj) - 5, 1):min(x_ind(jj) + 5, x_size),...
                max(y_ind(jj) - 5, 1):min(y_ind(jj) + 5, y_size), max(z_ind(jj) - 5, 1):min(z_ind(jj) + 5, z_size));
            if max(distance_win(:)) == distance_temp(x_ind(jj), y_ind(jj), z_ind(jj)) && max(distance_win(:)) > 0.5
                distance_max(x_ind(jj), y_ind(jj), z_ind(jj)) = distance_temp(x_ind(jj), y_ind(jj), z_ind(jj));
                % region growing under certain std
                max_std = std(distance_win(distance_win>0));
                thr = distance_temp(x_ind(jj), y_ind(jj), z_ind(jj)) - 2*max_std;
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
end

% x = distance_max(distance_max>0);
% fore_tip = zeros(size(fore_all));
% fore_tip(distance_max>0) = x - 0.5;
% 
% fore_tip = imdilate(fore_tip,ones(3,3,3));
% im = zeros(x_size, y_size, 3, z_size);
% for zz = 1:z_size
%     im(:,:,1,zz) = fore_temp(:,:,zz)*128;
%     im(:,:,2,zz) = fore_temp(:,:,zz)*128;
%     im(:,:,3,zz) = fore_temp(:,:,zz)*128;
% end
% fore_tip = ceil(fore_tip*255/max(fore_tip(:))-1e-8) + 1;
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

% x = distance_max(distance_max>0);
fore_tip = zeros(size(fore_all));
% fore_tip(distance_max>0) = x - 0.5;

com = bwconncomp(distance_max>0);
for ii = 1:com.NumObjects
    [x_ind, y_ind ,z_ind] = ind2sub(size(fore_all),com.PixelIdxList{ii});
    distance_temp = distance_max(com.PixelIdxList{ii});
    distance_norm = distance_temp/sum(distance_temp);
%     fore_tip(round(sum(x_ind.*distance_norm)),round(sum(y_ind.*distance_norm)),round(sum(z_ind.*distance_norm))) = mean(distance_temp)-0.5;
    fore_tip(com.PixelIdxList{ii}) = mean(distance_temp) - 0.5;
end
% fore_tip = imdilate(fore_tip,strel('sphere',4));
fore_tip = imdilate(fore_tip,ones(3,3,3));
im = zeros(x_size, y_size, 3, z_size);
for zz = 1:z_size
    im(:,:,1,zz) = fore_temp(:,:,zz)*128;
    im(:,:,2,zz) = fore_temp(:,:,zz)*128;
    im(:,:,3,zz) = fore_temp(:,:,zz)*128;
end
fore_tip = ceil(fore_tip*255/max(fore_tip(:))-1e-8) + 1;
color = colormap;
for ii = 1:x_size
    for jj = 1:y_size
        for kk = 1:z_size
            if fore_tip(ii,jj,kk) > 1
                im(ii,jj,1,kk) = color(fore_tip(ii,jj,kk),1)*255;
                im(ii,jj,2,kk) = color(fore_tip(ii,jj,kk),2)*255;
                im(ii,jj,3,kk) = color(fore_tip(ii,jj,kk),3)*255;
            end
        end
    end
end

tifwrite(uint8(im),'.\temp_data\tip3');

