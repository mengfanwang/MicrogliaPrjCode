clc;clear;close all;

dbstop if error 
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\redetection_info_v2';
% mkdir([file_name '\tip_distance']);
mkdir(fullfile(file_name, 'tipDetection'));

% neighbor = cat(3, [0 0 0; 0 1 0; 0 0 0],...
%                   [0 1 0; 1 1 1; 0 1 0],...
%                   [0 0 0; 0 1 0; 0 0 0]);
neighbor = ones(3,3);
step = 20;
             
for tt = 2:34
ind = num2str(tt);
% ind = ind(2:4);
% fore_all = load([file_name '\foreground_c2\' ind '.mat']);
% fore_all = load([file_name '\foreground_c2_reconnect_iter1\' ind '.mat']);
% fore_all = logical(fore_all.foreground);
fore_all = tifread(fullfile(file_name, 'foreground', [num2str(tt) '.tif']));
fore_all = logical(fore_all);



% remove small regions
com = bwconncomp(fore_all);
temp = zeros(size(fore_all));
com_size = cellfun(@length,com.PixelIdxList);
for ii = 1:com.NumObjects
    if com_size(ii) > 20 % small region threshold
        temp(com.PixelIdxList{ii}) = fore_all(com.PixelIdxList{ii});
    end
end
fore_all = temp;

fore_all = double(fore_all);
bound_temp = fore_all - imerode(fore_all, neighbor);

fore_all = padarray(fore_all,[step step step],'symmetric','both');
bound_all = zeros(size(fore_all));
bound_all(step+1:end-step,step+1:end-step,step+1:end-step) = bound_temp;
fore_ind = find(fore_all);
node_num = length(fore_ind);
fore_all(fore_ind) = 1:node_num;
bound_ind = find(bound_all);
bound_num = length(bound_ind);

tic;
distance_map = geodesic_v3(node_num, fore_all, bound_num, bound_all,1:step, [1 1 2.168], 8);
toc
distance_map = distance_map(step+1:end-step,step+1:end-step,step+1:end-step,:);
% save([file_name '\tip_distance\' ind '.mat'],'distance_map');
% mkdir([file_name '\tip_distance\' ind]);
save(fullfile(file_name, 'tipDetection', [ind '.mat']), 'distance_map');
mkdir(fullfile(file_name, 'tipDetection', ind));
for ii = 1:step
%     zz = num2str(1000+ii);
%     zz = zz(2:4);
%     writeColormap(distance_map(:,:,:,ii),[file_name '\tip_distance\' ind '\' num2str(ii)]);
    writeColormap(distance_map(:,:,:,ii),fullfile(file_name, 'tipDetection', ind, num2str(ii)));
end
end


function writeColormap(distance_map, filename)
distance_map(distance_map>10) = 10;
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

