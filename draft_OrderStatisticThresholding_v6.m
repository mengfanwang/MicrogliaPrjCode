clc;clear;close all;
dbstop if error;
% dbstop at 100 if ismember(sub2ind([512 512 61],440,421,1),com_element) 
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
mkdir([file_name '\foreground']);
load([file_name '\data_c2']);
load([file_name '\stabilizeFunction_c2']);
data_all = double(data2);
[x,y,z,t] = size(data_all);
std = sqrt(max(variance));


% This method will use bottom-up stragety. Every time we threshold 
% intensity from bottom to up, and select isolated
% highest-zscore regions. We repeat this step until no significant area.
% Smoothed data will be used.

% it will constraint from both min size and max size

for tt = 9:9
    tt
    data = data_all(:,:,:,tt);
    data = interp1(0:255,stabilizeFunction,data);
    
% 10 neighbor connecivty
% conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
%               [0 1 0; 1 1 1; 0 1 0],...
%               [0 0 0; 0 1 0; 0 0 0]);
conn = ones(3,3);
gaussian_sigma = 0.6;
z_thre = 2;
minsize = 3;
maxsize = inf;

data_smoothed = imgaussfilt3(data,gaussian_sigma);

% kern = cat(3, [0 0 0; 0 -1/6 0; 0 0 0],...
%               [0 -1/6 0; -1/6 2 -1/6; 0 -1/6 0],...
%               [0 0 0; 0 -1/6 0; 0 0 0]);
% data_smoothed = imfilter(data,kern);
% ind = num2str(1000+tt); 
% ind = [file_name '\data_sharped\' ind(2:4)];
% tifwrite(uint8(data_sharped), ind);


x_direction = zeros(1,27);
y_direction = zeros(1,27);
z_direction = zeros(1,27);
ind_dir = 1;
for xx = -1:1
    for yy = -1:1
        for zz = -1:1
            x_direction(ind_dir) = xx;
            y_direction(ind_dir) = yy;
            z_direction(ind_dir) = zz;
            ind_dir = ind_dir + 1;
        end
    end
end

components_cell = cell(25,2);
foreground = zeros(x,y,z);
tic;
for iter = 1:1000
    iter
    data_nan = data;
    data_nan(foreground==1) = nan;
    score_map = zeros(x,y,z);
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
                for dir_ind = 1:27
                    x_neighbor = min(max(x_position + x_direction(dir_ind),1),x); 
                    y_neighbor = min(max(y_position + y_direction(dir_ind),1),y); 
                    z_neighbor = min(max(z_position + z_direction(dir_ind),1),z); 
                    pix1 = sub2ind([x,y,z],x_neighbor,y_neighbor,z_neighbor);
                    pix1 = pix1(~ismember(pix1,com_element));
                    pix1 = pix1(~ismember(pix1,com_neighbor));
                    com_neighbor = [com_neighbor;pix1];
                end
                % second layer neighbood 
                com_neighbor2 = [];
                [x_position,y_position,z_position] = ind2sub([x,y,z],com_neighbor);
                for dir_ind = 1:27
                    x_neighbor = min(max(x_position + x_direction(dir_ind),1),x); 
                    y_neighbor = min(max(y_position + y_direction(dir_ind),1),y); 
                    z_neighbor = min(max(z_position + z_direction(dir_ind),1),z); 
                    pix1 = sub2ind([x,y,z],x_neighbor,y_neighbor,z_neighbor);
                    pix1 = pix1(~ismember(pix1,com_element));
                    pix1 = pix1(~ismember(pix1,com_neighbor));
                    pix1 = pix1(~ismember(pix1,com_neighbor2));
                    com_neighbor2 = [com_neighbor2;pix1];
                end

                M = numel(com_element);
                N = numel(com_neighbor2);
                L = mean(data_nan(com_element)) - mean(data_nan(com_neighbor2),'omitnan');

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
        end
        if iter == 1
            components_cell(thre,1) = {components.PixelIdxList};
            components_cell(thre,2) = {z_scores};
        end
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
% ind = [file_name '\foreground_original\' ind(2:4)];
ind = [file_name '\temp\' ind(2:4) '_' num2str(gaussian_sigma) '_ori'];
tifwrite(uint8(foreground*255), ind);
save([ind '.mat'],'foreground','components_cell');
% bfsave(uint8(max(foreground*255,data)),['.\BinaryData\before_vst_mask.tiff']);
end