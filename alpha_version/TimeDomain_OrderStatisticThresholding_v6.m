clc;clear;close all;
dbstop if error;
% dbstop at 100 if ismember(sub2ind([512 512 61],440,421,1),com_element) 
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
mkdir([file_name '\foreground_iter1']);
load([file_name '\data_c2']);
load([file_name '\stabilizeFunction_c2']);
data_all = double(data2);
[x,y,z,t] = size(data_all);
std = sqrt(max(variance)/3);

data_all = interp1(0:255,stabilizeFunction,data_all);
% preprocessing. data are mapped to the previous and last frames to get a
% better estimation
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
% for tt = 1:t
%     ind = num2str(1000+tt);
%     ind = ind(2:4);
%     tifwrite(uint8(data_all(:,:,:,tt)), ['D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\time_domain_segmentation\' ind]);
% end


% This method will use bottom-up stragety. Every time we threshold 
% intensity from bottom to up, and select isolated
% highest-zscore regions. We repeat this step until no significant area.
% Smoothed data will be used.

% it will constraint from both min size and max size


for tt = 1:t
    tt
    data = data_all(:,:,:,tt);
    z_map = nan(size(data));
%     data = interp1(0:255,stabilizeFunction,data);
    
% 10 neighbor connecivty
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);
% conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
%               [1 1 1; 1 1 1; 1 1 1],...
%               [0 0 0; 0 1 0; 0 0 0]);
% conn = cat(3, [0 1 0; 1 1 1; 0 1 0],...
%               [1 1 1; 1 1 1; 1 1 1],...
%               [0 1 0; 1 1 1; 0 1 0]);
% conn = ones(3,3);
gaussian_sigma = 0.8;
z_thre = 2;
minsize = 3;
maxsize = 10000;

data_smoothed = imgaussfilt3(data,gaussian_sigma);

% kern = cat(3, [0 0 0; 0 -1/6 0; 0 0 0],...
%               [0 -1/6 0; -1/6 2 -1/6; 0 -1/6 0],...
%               [0 0 0; 0 -1/6 0; 0 0 0]);
% data_smoothed = imfilter(data,kern);
% ind = num2str(1000+tt); 
% ind = [file_name '\data_sharped\' ind(2:4)];
% tifwrite(uint8(data_sharped), ind);


% x_direction = zeros(1,27);
% y_direction = zeros(1,27);
% z_direction = zeros(1,27);
% ind_dir = 1;
% for xx = -1:1
%     for yy = -1:1
%         for zz = -1:1
%             x_direction(ind_dir) = xx;
%             y_direction(ind_dir) = yy;
%             z_direction(ind_dir) = zz;
%             ind_dir = ind_dir + 1;
%         end
%     end
% end
x_direction = [1 -1 0 0 0 0];
y_direction = [0 0 1 -1 0 0];
z_direction = [0 0 0 0 1 -1];

components_cell = cell(100,1);
foreground = zeros(x,y,z);
tic;
for iter = 1:1000
    iter
    data_nan = data;
    data_nan(foreground==1) = nan;
    score_map = zeros(x,y,z);
    z_map_temp = nan(x,y,z);
    components_cell{iter} = cell(25,2);
    for thre = 1:25
        binaryData = zeros(x,y,z);
        binaryData(data_smoothed<= thre*10) = 0;
        binaryData(data_smoothed>  thre*10) = 1;
%         binaryData = imopen(binaryData, conn);
        
        binaryData(foreground == 1) = 0;
%         binaryData = binaryData.*open_data; % open operation

        components = bwconncomp(binaryData,conn);
        components_len = length(components.PixelIdxList);
        z_scores = zeros(components_len,1);
        for com_ind = 1:components_len
    %        com_ind    
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
                
%                 % second layer neighbood 
%                 com_neighbor2 = [];
%                 [x_position,y_position,z_position] = ind2sub([x,y,z],com_neighbor);
%                 for dir_ind = 1:length(x_direction)
%                     x_neighbor = min(max(x_position + x_direction(dir_ind),1),x); 
%                     y_neighbor = min(max(y_position + y_direction(dir_ind),1),y); 
%                     z_neighbor = min(max(z_position + z_direction(dir_ind),1),z); 
%                     pix1 = sub2ind([x,y,z],x_neighbor,y_neighbor,z_neighbor);
%                     pix1 = pix1(~ismember(pix1,com_element));
%                     pix1 = pix1(~ismember(pix1,com_neighbor));
%                     pix1 = pix1(~ismember(pix1,com_neighbor2));
%                     com_neighbor2 = [com_neighbor2;pix1];
%                 end
%                 com_neighbor2 = unique(com_neighbor2);
%                 
%                 % first layer element
%                 com_element2 = [];
%                 [x_position,y_position,z_position] = ind2sub([x,y,z],com_neighbor);
%                 for dir_ind = 1:length(x_direction)
%                     x_neighbor = min(max(x_position + x_direction(dir_ind),1),x); 
%                     y_neighbor = min(max(y_position + y_direction(dir_ind),1),y); 
%                     z_neighbor = min(max(z_position + z_direction(dir_ind),1),z); 
%                     pix1 = sub2ind([x,y,z],x_neighbor,y_neighbor,z_neighbor);
%                     pix1 = pix1(ismember(pix1,com_element));
%                     com_element2 = [com_element2;pix1];
%                 end
%                 com_element2 = unique(com_element2);
%                 com_element = com_element2;
                
                
                L = mean(data_nan(com_element)) - mean(data_nan(com_neighbor));
                largeGroup = data_nan(com_element);
                smallGroup = data_nan(com_neighbor);
                [mu, sigma] = ordStatApproxKsec(largeGroup, smallGroup);
                z_scores(com_ind) = (L/std - mu)/sigma;
            end
        end
        for com_ind = 1:components_len 
            if z_scores(com_ind)>z_thre
                score_map(components.PixelIdxList{com_ind}) = max(z_scores(com_ind),score_map(components.PixelIdxList{com_ind}));
            end
            z_map_temp(components.PixelIdxList{com_ind}) = max(z_scores(com_ind),z_map_temp(components.PixelIdxList{com_ind}));
        end
        
        components_cell{iter}(thre,1) = {components.PixelIdxList};
        components_cell{iter}(thre,2) = {z_scores};
    end
    if sum(score_map(:)) == 0
        break
    end
    foreground = foreground + imregionalmax(score_map,26);
    z_map(imregionalmax(score_map,26)) = score_map(imregionalmax(score_map,26));
    tifwrite(uint8(z_map*5), [file_name '\temp\z_map\' num2str(iter)]);
end
toc
z_map(isnan(z_map)) = z_map_temp(isnan(z_map));

%% remove tiny objects
components = bwconncomp(foreground,conn);
foreground = zeros(size(foreground));
for ii = 1:components.NumObjects
    if length(components.PixelIdxList{ii}) > 20
        foreground(components.PixelIdxList{ii}) = 1;
    end
end

% save(['.\BinaryData\v6\z_score_n2\' num2str(thre) '.mat'],'components_cell');
ind = num2str(1000+tt); 
ind = [file_name '\foreground_iter1_simulation\' ind(2:4)];
tifwrite(uint8(foreground*255), ind);
% tifwrite(uint8(foreground*255), [ind '_' num2str(z_thre) '_' num2str(sum(conn(:))-1) '_' num2str(gaussian_sigma) '_1layer']);
components_cell = components_cell(1:iter);
save([ind '.mat'],'foreground','components_cell','z_map');
% save([file_name '\temp\z_map\' ind '.mat'],'foreground','z_map');
end

