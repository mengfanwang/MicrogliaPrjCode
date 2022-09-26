clc;clear;close all;

soma_ind = 4;
file_name = '091919 slice1 cortex 4mM isoflurane';

mkdir([file_name '\cell']);
mkdir([file_name '\result']);
load([file_name '\foreground_reg']);
load([file_name '\soma_reg']);


[x_min, y_min, z_min] = deal(inf);
[x_max, y_max, z_max] = deal(-inf);

% find min/max region
for tt = 1:size(foreground_reg,4)
    soma = soma_reg(:,:,:,tt);
    foreground = foreground_reg(:,:,:,tt);
    [x_size, y_size, z_size] = size(foreground);

    % find soma
    soma_com = bwconncomp(soma,26);
    if tt == 1
        soma = zeros(size(soma));
        soma(soma_com.PixelIdxList{soma_ind}) = 1;
        soma = logical(soma);
        soma_prev = find(soma);
    else
        for ii = 1:length(soma_com.PixelIdxList)
            if any(ismember(soma_com.PixelIdxList{ii},soma_prev))
                soma = zeros(size(soma));
                soma(soma_com.PixelIdxList{ii}) = 1;
                soma = logical(soma);
                soma_prev = find(soma);
                soma_label = 1;
            end
        end
    end
    
    % find foreground
    foreground_bg = foreground;
    foreground_com = bwconncomp(foreground,26);
    for ii = 1:length(foreground_com.PixelIdxList)
        if any(ismember(foreground_com.PixelIdxList{ii},soma_prev))
            foreground = zeros(size(foreground));
            foreground(foreground_com.PixelIdxList{ii}) = 1;
            foreground = logical(foreground);
        end
    end
    
    
    [x_ind, y_ind, z_ind] = ind2sub(size(foreground),find(foreground));
    x_min = min(x_min,min(x_ind));
    y_min = min(y_min,min(y_ind));
    z_min = min(z_min,min(z_ind));
    x_max = max(x_max,max(x_ind));
    y_max = max(y_max,max(y_ind));
    z_max = max(z_max,max(z_ind));
end

for tt = 1:size(foreground_reg,4)
    ind = num2str(1000+tt);
    ind = ind(2:4);
    foreground = foreground_reg(x_min:x_max,y_min:y_max,z_min:z_max);
    tifwrite(uint8(foreground*255),[file_name '\cell\' ind]);
end

for tt = 1:size(foreground_reg,4)
    tt
    tic;
    ind = num2str(1000+tt);
    ind = ind(2:4);
    soma = soma_reg(:,:,:,tt);
    foreground = foreground_reg(:,:,:,tt);

    % find soma
    soma_com = bwconncomp(soma,26);
    if tt == 1
        soma = zeros(size(soma));
        soma(soma_com.PixelIdxList{soma_ind}) = 1;
        soma = logical(soma);
        soma_prev = find(soma);
    else
        for ii = 1:length(soma_com.PixelIdxList)
            if any(ismember(soma_com.PixelIdxList{ii},soma_prev))
                soma = zeros(size(soma));
                soma(soma_com.PixelIdxList{ii}) = 1;
                soma = logical(soma);
                soma_prev = find(soma);
            end
        end
    end
    
    % find foreground
    foreground_bg = foreground;
    foreground_com = bwconncomp(foreground,26);
    for ii = 1:length(foreground_com.PixelIdxList)
        if any(ismember(foreground_com.PixelIdxList{ii},soma_prev))
            foreground = zeros(size(foreground));
            foreground(foreground_com.PixelIdxList{ii}) = 1;
            foreground = logical(foreground);
        end
    end
    foreground_bg = foreground_bg - foreground;
    
   
    soma = soma(x_min:x_max,y_min:y_max,z_min:z_max);
    foreground = foreground(x_min:x_max,y_min:y_max,z_min:z_max);
    foreground_bg = foreground_bg(x_min:x_max,y_min:y_max,z_min:z_max);
    [result{tt}, ~, ~, ~, Node{tt}, G{tt}] = treeBuilding(foreground,soma);
    cell{tt} = foreground;
    cell_bg{tt} = foreground + foreground_bg;
    
    [x_size, y_size, z_size] = size(foreground_bg);
    for zz = 1:z_size
        result{tt}(:,:,1,zz) = result{tt}(:,:,1,zz) + foreground_bg(:,:,zz)*128;
        result{tt}(:,:,2,zz) = result{tt}(:,:,2,zz) + foreground_bg(:,:,zz)*128;
        result{tt}(:,:,3,zz) = result{tt}(:,:,3,zz) + foreground_bg(:,:,zz)*128;
    end
    toc
    tifwrite(uint8(result{tt}),[file_name '\result\' ind]);
end
save([file_name '\result'],'result','-v7.3');
save([file_name '\graph'],'Node','G');
save([file_name '\cell'],'cell','-v7.3');
save([file_name '\cell_bg'],'cell_bg','-v7.3');

