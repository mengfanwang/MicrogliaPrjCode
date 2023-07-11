clc;clear;close all;
x_size = 512; y_size = 511; z_size = 61;
% the last column is deleted becasue of the data property
% in the original data, the last column is almost empty
dbstop if error 
for tt = 1:1
    tt
ind = num2str(1000+tt);
ind = ind(2:4);

win_size = 15;
neighbor = cat(3, [0 0 0; 0 1 0; 0 0 0],...
                  [0 1 0; 1 1 1; 0 1 0],...
                  [0 0 0; 0 1 0; 0 0 0]);
% calculate distance
step = (win_size-1)/2;
center = [(win_size+1)/2 (win_size+1)/2 (win_size+1)/2];
center_idx = sub2ind([win_size win_size win_size],center(1),center(2),center(3));

% load whole data
fore_all = load(['D:\dropbox\foreground\' ind '.mat']);
fore_all = logical(fore_all.labelMap);
fore_all = fore_all(:,1:y_size,:);
com = bwconncomp(fore_all);
distance_map = nan(x_size,y_size,z_size);
distance = nan;
tic;
for ii = 1:length(com.PixelIdxList)
    
    fore_temp = zeros(size(fore_all));
    fore_temp(com.PixelIdxList{ii}) = 1;
    fore_temp = logical(fore_temp);
    fore_bound = fore_temp - imerode(fore_temp, neighbor);
    [x_ind,y_ind,z_ind] = ind2sub(size(fore_bound),find(fore_bound));
    edge_list = [x_ind y_ind z_ind];

    for jj = 1:size(edge_list,1)
        if edge_list(jj,1) < (win_size+1)/2 || edge_list(jj,2) < (win_size+1)/2 || edge_list(jj,3) < (win_size+1)/2 ||...
           edge_list(jj,1) >  x_size - (win_size-1)/2 || edge_list(jj,2) >  y_size - (win_size-1)/2 || edge_list(jj,3) >  z_size - (win_size-1)/2
            continue
        end
        win = fore_temp(edge_list(jj,1)-step:edge_list(jj,1)+step,edge_list(jj,2)-step:edge_list(jj,2)+step,edge_list(jj,3)-step:edge_list(jj,3)+step);
        dist_map = bwdistgeodesic(win,center_idx,'quasi-euclidean');
        [x_ind,y_ind, z_ind] = ind2sub(size(dist_map),find(dist_map > step-0.5 & dist_map < step+0.5));
        line_ind =[x_ind y_ind z_ind];
        if isempty(x_ind)
            distance = nan;
        elseif length(x_ind) == 1
            distance = norm(center - line_ind);
        elseif length(x_ind) == 2
            distance = calculatePoint2LineDist(center,line_ind(1,:),line_ind(2,:));
        else

            try
                [conv_ind, v1] = convhull(x_ind,y_ind,z_ind,'Simplify',true);
                conv_ind = unique(conv_ind);
                [x,distance] = findDistance(center',line_ind(conv_ind,:)');                    
%                     distance = findDistance_draft(center', line_ind');
            catch 
                [x,distance] = findDistance(center',line_ind');
            end
%                 distance2 = openGJK(center', line_ind');

        end
%               distance_map(edge_list(jj,1),edge_list(jj,2),edge_list(jj,3)) = distance;
%             if ~isempty(x_ind)
% %             distance = openGJK(center', line_ind');
%                 distance2 = openGJK_v2(line_ind'-center');
% %                 abs(distance - distance2) > 1e-2;
% % dlmwrite('C:\Users\mengf\Desktop\openGJK_v2\userP.dat',size(line_ind,1));
% % dlmwrite('C:\Users\mengf\Desktop\openGJK_v2\userP.dat',line_ind-center,'delimiter',' ','-append');
% 
%              distance_map(edge_list(jj,1),edge_list(jj,2),edge_list(jj,3)) = distance2;
%             end
    end
end
running_time(tt) = toc;
% save(['.\test_data\distance_' name '_' ind], 'distance_map');
% map_3 = max(distance_map,[],3);
% % imagesc(map_3);colorbar;
% imagesc(map_3,'AlphaData',~isnan(map_3));
% set(gca,'color',0*[1 1 1]);
% colorbar;
end
save(['.\test_data\time_' name],'running_time');

function [x, distance] = findDistance(p,A)
    num = size(A,2);
    H = 2*A'*A;
    f = -2*A'*p;
    Aeq = ones(1,num);
    beq = 1;
    lb = zeros(num,1);
    ub = ones(num,1);
    options = optimoptions('quadprog','Display','off');
    [x,fval] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
    distance = sqrt(abs(fval + p'*p));
end

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