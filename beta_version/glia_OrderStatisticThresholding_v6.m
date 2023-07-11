% clc;clear;close all;
file_name = 'D:\dropbox\Modify Series Data-4nd\SL4\Copy of 042519 slice4 cerebellum before_Modify Series';

% function glia_OrderStatisticThresholding_v6(file_name)
channel = 'c2';
load([file_name '\registration_' channel]);
data_all = data2;

mkdir([file_name '\foreground_' channel]);
load([file_name '\stabilizeFunction_' channel]);
[x,y,z,t] = size(data_all);
std = sqrt(max(variance));


% parameter
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [0 1 0; 1 1 1; 0 1 0],...
              [0 0 0; 0 1 0; 0 0 0]);
% conn = ones(3,3);
gaussian_sigma = 0.8;
z_thre = 2;
minsize = 3;
maxsize = 20000;
x_direction = [1 -1 0 0 0 0];
y_direction = [0 0 1 -1 0 0];
z_direction = [0 0 0 0 1 -1];
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

% This method will use bottom-up stragety. Every time we threshold 
% intensity from bottom to up, and select isolated
% highest-zscore regions. We repeat this step until no significant area.
% Smoothed data will be used.

% it will constraint from both min size and max size

for tt = 1:t
    tt
    data = data_all(:,:,:,tt);
    z_map = nan(size(data));
    data = interp1(0:255,stabilizeFunction,data);
%     data = data*1.2;
    data(data>255) = 255;


data_smoothed = imgaussfilt3(data,gaussian_sigma);


foreground = zeros(x,y,z);
tic;
for iter = 1:1000
%     iter
    data_nan = data;
    data_nan(foreground==1) = nan;
    score_map = zeros(x,y,z);
    z_map_temp = nan(x,y,z);
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
                
%                 L = mean(data_nan(com_element)) - mean(data_nan(com_neighbor));
%                 largeGroup = data_nan(com_element);
%                 smallGroup = data_nan(com_neighbor);
%                 [mu, sigma] = ordStatApproxKsec(largeGroup, smallGroup);
%                 z_scores(com_ind) = (L/std - mu)/sigma;
                
                M = numel(com_element);
                N = numel(com_neighbor);
                L = mean(data_nan(com_element)) - mean(data_nan(com_neighbor),'omitnan');
                pdfz = normpdf(norminv(N/(M+N),0,1));
                cdfz = normcdf(norminv(N/(M+N),0,1));
                mu = pdfz/(1-cdfz) + pdfz/cdfz;
                sigma = sqrt((1+norminv(N/(M+N),0,1)*pdfz/(1-cdfz)-(pdfz/(1-cdfz))^2)/M + (1-norminv(N/(M+N),0,1)*pdfz/cdfz-(pdfz/cdfz)^2)/N);
                z_scores(com_ind) = (L/std - mu)/sigma;
                
            end
        end
        for com_ind = 1:components_len 
            if z_scores(com_ind)>z_thre
                score_map(components.PixelIdxList{com_ind}) = max(z_scores(com_ind),score_map(components.PixelIdxList{com_ind}));
            end
            z_map_temp(components.PixelIdxList{com_ind}) = max(z_scores(com_ind),z_map_temp(components.PixelIdxList{com_ind}));
        end
        
    end
    if sum(score_map(:)) == 0
        break
    end
    foreground = foreground + imregionalmax(score_map,26);
%     bfsave(uint8(foreground*255),['.\BinaryData\v6\bottomup_05_' num2str(z_thre) '_' num2str(minsize) '_n2_0change1_255\' num2str(iter) '.tiff']);
end
toc
%% remove tiny objects
components = bwconncomp(foreground,conn);
foreground = zeros(size(foreground));
for ii = 1:components.NumObjects
    if length(components.PixelIdxList{ii}) > 100
        foreground(components.PixelIdxList{ii}) = 1;
    end
end

% save(['.\BinaryData\v6\z_score_n2\' num2str(thre) '.mat'],'components_cell');
ind = num2str(1000+tt); 
ind = [file_name '\foreground_' channel '\' ind(2:4)];
tifwrite(uint8(foreground*255), ind);
save([ind '.mat'],'foreground');
% bfsave(uint8(max(foreground*255,data)),['.\BinaryData\before_vst_mask.tiff']);
end