clc;clear;close all;
dbstop if error;
% dbstop at 100 if ismember(sub2ind([512 512 61],440,421,1),com_element) 
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
mkdir([file_name '\foreground_c1']);
load([file_name '\data_c1_processed']);
load([file_name '\stabilizeFunction_c1']);
data_all = double(data1);
[x,y,z,t] = size(data_all);
std = sqrt(max(variance));

% This method will use bottom-up stragety. Every time we threshold 
% intensity from bottom to up, and select isolated
% highest-zscore regions. We repeat this step until no significant area.
% Smoothed data will be used.

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
%     bfsave(uint8(foreground*255),['.\BinaryData\v6\bottomup_05_' num2str(z_thre) '_' num2str(minsize) '_n2_0change1_255\' num2str(iter) '.tiff']);
end
toc

% save(['.\BinaryData\v6\z_score_n2\' num2str(thre) '.mat'],'components_cell');
ind = num2str(1000+tt); 
ind = [file_name '\foreground_c1\' ind(2:4)];
tifwrite(uint8(foreground*255), ind);
save([ind '.mat'],'foreground');
% bfsave(uint8(max(foreground*255,data)),['.\BinaryData\before_vst_mask.tiff']);
end