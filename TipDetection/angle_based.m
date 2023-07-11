clc;clear;close all;
x_size = 512; y_size = 511; z_size = 61;
% the last column is deleted becasue of the data property
% in the original data, the last column is almost empty

win_size = 15;
[check_table, win_table] = checkVectorAngle(win_size);
for tt = 1:1
    tt
% tt = 1;
ind = num2str(1000+tt);
ind = ind(2:4);

% load whole data
fore_all = load(['D:\dropbox\foreground\' ind '.mat']);
fore_all = logical(fore_all.labelMap);
fore_all = fore_all(:,1:y_size,:);
com = bwconncomp(fore_all);
angle_map = nan(x_size,y_size,length(com.PixelIdxList));
for ii = 1:length(com.PixelIdxList)
    fore_temp = zeros(size(fore_all));
    fore_temp(com.PixelIdxList{ii}) = 1;
    fore_temp = max(fore_temp,[],3);
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

        % count histogram
        step = (win_size-1)/2;
        center = [(win_size+1)/2 (win_size+1)/2];
        
        for jj = 1:size(edge_list,1)
            if edge_list(jj,1) < (win_size+1)/2 || edge_list(jj,2) < (win_size+1)/2 ||...
                    edge_list(jj,1) >  x_size - (win_size-1)/2 || edge_list(jj,2) >  y_size - (win_size-1)/2
                continue
            end
            if jj <= step
                jj_previous = size(edge_list,1) - (step - jj);
                jj_after = jj + step;
                bound_check = sum(edge_list(jj_previous:end,3)) + sum(edge_list(1:jj_after,3));
            elseif jj > size(edge_list,1) - step
                jj_previous = jj - step;
                jj_after = step - (size(edge_list,1) - jj);
                bound_check = sum(edge_list(jj_previous:end,3)) + sum(edge_list(1:jj_after,3));
            else
                jj_previous = jj - step;
                jj_after = jj + step;
                bound_check = sum(edge_list(jj_previous:jj_after,3));
            end
            if bound_check == 0
                p1 = edge_list(jj_previous,1:2) - edge_list(jj,1:2) + center;
                p2 = edge_list(jj_after,1:2) - edge_list(jj,1:2) + center;
                if all(p1 == center) || all(p2 == center)
                    continue
                end
                info = check_table{p1(1),p1(2)}(win_table(p2(1),p2(2)));
                table1 = info.table1;
                table2 = info.table2;
                angle = info.angle;
                win = fore_temp(edge_list(jj,1)-step:edge_list(jj,1)+step,edge_list(jj,2)-step:edge_list(jj,2)+step);
                if info.angle < 1 - 1e-8 && info.angle > 1e-8 && ...
                        sum(table1.*win/sum(table1(:)),'all') < sum(table2.*win/sum(table2(:)),'all')
                   angle = 2- angle;
                end
                angle_map(edge_list(jj,1),edge_list(jj,2),ii) = angle;
            end 
        end
    end
    
end
data{tt} = angle_map(~isnan(angle_map));
end
figure(1);
histogram(angle_map);
figure(2);
map_3 = max(angle_map,[],3);
imagesc(map_3,'AlphaData',~isnan(map_3));
set(gca,'color',0*[1 1 1]);
oldcmap = colormap;
colormap( flipud(oldcmap) );
colorbar;


function [check_table, win_table] = checkVectorAngle(win_size)
% comapre which part is foreground based on the fore area percentage
center = [(win_size+1)/2 (win_size+1)/2];

tic;
check_table = cell(win_size,win_size);
win_table = zeros(win_size, win_size);
ind = 0;
for ii = 1:win_size
    for jj = 1:win_size
        if all([ii jj] == center)
            continue
        end
        ind = ind + 1;
        win_table(ii,jj) = ind;
    end
end
for mm = 1:win_size
    for nn = 1:win_size
        if all([mm nn] == center)
            continue
        end
        ind = 0;
        for kk = 1:win_size
            for ll = 1:win_size
                if all([kk ll] == center)
                    continue
                end
                p1 = [mm nn];
                p2 = [kk ll];
                v1 = p1 - center;
                v2 = p2 - center;
                angle = real(acos((v1*v2')/norm(v1)/norm(v2)))/pi;    % angle/pi
                n1 = [v1(2) -v1(1)];
                n2 = [-v2(2) v2(1)];                            % in case of angle = -pi
                if n1*v2'< 0 && angle ~= 1
                    n1 = -n1;
                end
                if n2*v1'< 0 && angle ~= 1
                    n2 = -n2;
                end
                table1 = zeros(win_size,win_size);
                table2 = zeros(win_size,win_size);
                for ii = 1:win_size
                    for jj = 1:win_size
                        v_test = [ii-center(1) jj-center(2)];
                        if v_test*v1' == norm(v_test)*norm(v1) || v_test*v2' == norm(v_test)*norm(v2)
                            continue
                        end
                        if all(v_test == [0 0])
                            continue
                        end
                        if v_test*n1' > 0 && v_test*n2'>0
                            table1(ii,jj) = 1;
                        end
                        if v_test*n1' < 0 || v_test*n2'<0
                            table2(ii,jj) = 1;
                        end
                    end
                end
                ind = ind + 1;
                check_table{mm,nn}(ind).p2 = p2;
                check_table{mm,nn}(ind).angle = angle;
                check_table{mm,nn}(ind).table1 = table1;
                check_table{mm,nn}(ind).table2 = table2;
            end
        end
    end
end
toc
end

