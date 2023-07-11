clc;clear;close all;
x_size = 512; y_size = 511; z_size = 61;
% the last column is deleted becasue of the data property
% in the original data, the last column is almost empty

tt = 1;
ind = num2str(1000+tt);
ind = ind(2:4);

win_size = 15;

% load whole data
fore_all = load(['D:\dropbox\foreground\' ind '.mat']);
fore_all = logical(fore_all.labelMap);
fore_all = fore_all(:,1:y_size,:);
com = bwconncomp(fore_all);
distance_map = nan(x_size,y_size,length(com.PixelIdxList));
tic;
for ii = 1:length(com.PixelIdxList)
    fore_temp = zeros(size(fore_all));
    fore_temp(com.PixelIdxList{ii}) = 1;
    fore_temp = max(fore_temp,[],3);
    fore_temp = logical(fore_temp);
    if sum(fore_temp(:)) > 100
        edge_list = bwboundaries(fore_temp,8,'noholes');
        edge_list = edge_list{1}(1:end-1,:);
        
        % find boundary pixels
        bound_pix = zeros(size(edge_list,1),1);
        for jj = 1:size(edge_list,1)
            if jj == 1
                jj_previous = size(edge_list,1);
                jj_after = jj + 1;
            elseif jj == size(edge_list,1)
                jj_previous = jj - 1;
                jj_after = 1;
            else
                jj_previous = jj - 1;
                jj_after = jj + 1;
            end
            
            if edge_list(jj,1) == 1 && edge_list(jj_previous,1) == 1 && edge_list(jj_after,1) == 1
                bound_pix(jj) = 1;
            end
            if edge_list(jj,1) == x_size && edge_list(jj_previous,1) == x_size && edge_list(jj_after,1) == x_size
                bound_pix(jj) = 1;
            end
            if edge_list(jj,2) == 1 && edge_list(jj_previous,2) == 1 && edge_list(jj_after,2) == 1
                bound_pix(jj) = 1;
            end
            if edge_list(jj,2) == y_size && edge_list(jj_previous,2) == y_size && edge_list(jj_after,2) == y_size
                bound_pix(jj) = 1;
            end
        end
        
        edge_list = [edge_list bound_pix];
        edge_im = zeros(x_size,y_size,3);
        for jj = 1:size(edge_list,1)
            edge_im(edge_list(jj,1),edge_list(jj,2),:) = edge_im(edge_list(jj,1),edge_list(jj,2)) + 127;
            if edge_list(jj,3)
                edge_im(edge_list(jj,1),edge_list(jj,2),2:3) = 0;
            end
        end
%         imshow(uint8(edge_im));
%         imwrite(uint8(edge_im),['.\temp\edge\' num2str(ii) '.tif']);

        % calculate distance
        step = (win_size-1)/2;
        center = [(win_size+1)/2 (win_size+1)/2];
        
        for jj = 1:size(edge_list,1)
            if edge_list(jj,1) < (win_size+1)/2 || edge_list(jj,2) < (win_size+1)/2 ||...
                    edge_list(jj,1) >  x_size - (win_size-1)/2 || edge_list(jj,2) >  y_size - (win_size-1)/2
                continue
            end
            win = fore_temp(edge_list(jj,1)-step:edge_list(jj,1)+step,edge_list(jj,2)-step:edge_list(jj,2)+step);
            dist_map = bwdistgeodesic(win,center(1),center(2),'quasi-euclidean');
            [x_ind,y_ind] = find(dist_map > step-1 & dist_map < step+1);
            line_ind =[x_ind y_ind];
            if isempty(x_ind)
                min_dist = nan;
            elseif length(x_ind) == 1
                min_dist = norm(center - line_ind);
            else
                min_dist = inf;
                combination = nchoosek(1:length(x_ind),2);
                for kk = 1:size(combination,1)
                    min_dist = min(min_dist, calculatePoint2LineDist(center,line_ind(combination(kk,1),:),line_ind(combination(kk,2),:)));
                end
            end
            distance_map(edge_list(jj,1),edge_list(jj,2),ii) = length(x_ind);%min_dist;
        end
    end
end
toc
map_3 = max(distance_map,[],3);
% imagesc(map_3);colorbar;
imagesc(map_3,'AlphaData',~isnan(map_3));
set(gca,'color',0*[1 1 1]);
colorbar;

function distance = calculatePoint2LineDist(p,p1,p2)
    p1_p = p - p1;
    p1_p2 = p2 - p1;
    theta = sum(p1_p.*p1_p2)/sum(p1_p2.^2);
    if theta < 0
        pp = p1;
    elseif theta > 1
        pp = p2;
    else
        pp = p1 + theta*p1_p2;
    end
    distance = norm(p - pp);
end
