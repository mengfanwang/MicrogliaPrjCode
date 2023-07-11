clc;clear;close all;
dbstop if error

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';
mkdir([file_name '\tip_tracking']);
load([file_name 'data_c2.mat']);
data_all = data2;
load([file_name 'tip_info_reconnect.mat']);
[x, y, z, t] = size(data_all);
addpath ./CINDA


%%% parameter setting %%%
v_max = 20;
d_max = 15;


% multiframe graph
num_tip = cellfun(@length,xCoord);
% cum_tip = cumsum(num_tip);
detection_arcs = zeros(sum(num_tip),4);
node_table = zeros(sum(num_tip),4);
node_matrix = zeros(size(data_all));
transition_arcs = zeros(sum(num_tip(1:t-1).*num_tip(2:t)),3);
transition_ind = 0;
for tt = 1:t
    tt
    for jj = 1:num_tip(tt)
        ind = sum(num_tip(1:tt-1))+jj;
        detection_arcs(ind,:) = [ind v_max/2 v_max/2 -v_max];
        node_table(ind,:) = [tt xCoord{tt}(jj) yCoord{tt}(jj) zCoord{tt}(jj)];
        node_matrix(xCoord{tt}(jj),yCoord{tt}(jj),zCoord{tt}(jj),tt) = ind;
        if tt ~= t
            for kk = 1:num_tip(tt+1)
                dist_temp = sqrt((xCoord{tt}(jj)-xCoord{tt+1}(kk))^2 + (yCoord{tt}(jj)-yCoord{tt+1}(kk))^2 + (zCoord{tt}(jj)-zCoord{tt+1}(kk))^2);
                if dist_temp < d_max
                    if xCoord{tt}(jj) >= 8 && yCoord{tt}(jj) >= 8 && zCoord{tt}(jj)>=8 && ...
                        xCoord{tt}(jj) <= x-7 && yCoord{tt}(jj) <= y-7 && zCoord{tt}(jj) <= z-7 && ...
                        xCoord{tt+1}(kk) >= 8 && yCoord{tt+1}(kk) >= 8 && zCoord{tt+1}(kk)>=8 && ...
                        xCoord{tt+1}(kk) <= x-7 && yCoord{tt+1}(kk) <= y-7 && zCoord{tt+1}(kk) <= z-7 
                    
                        transition_ind = transition_ind + 1;
                        local1 = data_all(xCoord{tt}(jj)-7:xCoord{tt}(jj)+7, yCoord{tt}(jj)-7:yCoord{tt}(jj)+7, zCoord{tt}(jj)-7:zCoord{tt}(jj)+7,tt);
                        local2 = data_all(xCoord{tt+1}(kk)-7:xCoord{tt+1}(kk)+7, yCoord{tt+1}(kk)-7:yCoord{tt+1}(kk)+7, zCoord{tt+1}(kk)-7:zCoord{tt+1}(kk)+7,tt+1);
                        local1_pattern = zeros(5,5,5);
                        local2_pattern = zeros(5,5,5);
                        for xx = 1:5
                            for yy = 1:5
                                for zz = 1:5
                                    local1_pattern(xx,yy,zz) = mean(local1(xx*3-2:xx*3,yy*3-2:yy*3,zz*3-2:zz*3),'all');
                                    local2_pattern(xx,yy,zz) = mean(local2(xx*3-2:xx*3,yy*3-2:yy*3,zz*3-2:zz*3),'all');
                                end
                            end
                        end
                        diff_temp = sum(abs(local1_pattern-local2_pattern)/255, 'all');
                        transition_arcs(transition_ind,:) = [ind sum(num_tip(1:tt))+kk diff_temp];
                    end
                end
            end
        end
    end
end
transition_arcs = transition_arcs(1:transition_ind,:);

[trajectories, costs] = mcc4mot(detection_arcs,transition_arcs);
trajectories(cellfun(@length,trajectories)<5) = [];

% draw 3D line
im = zeros(x, y, 3, z, t);
fore_tip = zeros(x,y,z,t);
for ii = 1:length(trajectories)
    for jj = 1:length(trajectories{ii})
        tip1 = trajectories{ii}(jj);
        fore_tip(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),node_table(tip1,1)) = 1;
    end
end
tic;
for ii = 1:length(trajectories)
    ii
    for jj = 1:length(trajectories{ii})-1
        pix1 = trajectories{ii}(jj);
        pix2 = trajectories{ii}(jj+1);
        for tt = node_table(pix2,1):node_table(trajectories{ii}(end),1)
            headYXZ = [node_table(pix1,2) node_table(pix1,3) node_table(pix1,4)];
            tailYXZ = [node_table(pix2,2) node_table(pix2,3) node_table(pix2,4)];
            im(:,:,:,:,tt) = draw3Dline(im(:,:,:,:,tt), headYXZ, tailYXZ, 1, [0 255 0]);
        end
    end
end
toc
for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    foreground = load([file_name 'foreground_c2\' ind '.mat']);
    foreground = foreground.foreground;
    
    tip_dilate = imdilate(fore_tip(:,:,:,tt),strel('sphere',3));
    for zz = 1:z
        im(:,:,1,zz,tt) = im(:,:,1,zz,tt) + foreground(:,:,zz)*128;
        im(:,:,2,zz,tt) = im(:,:,2,zz,tt) + foreground(:,:,zz)*128 + tip_dilate(:,:,zz)*255;
        im(:,:,3,zz,tt) = im(:,:,3,zz,tt) + foreground(:,:,zz)*128 + tip_dilate(:,:,zz)*255;
    end
    
    tifwrite(uint8(im(:,:,:,:,tt)), [file_name '\tip_tracking_LBP\' ind ]);
end




a = 1;


for zz = 1:z
    im(:,:,1,zz) = foreground(:,:,zz)*128;
    im(:,:,2,zz) = foreground(:,:,zz)*128 + tip_dilate(:,:,zz)*255;
    im(:,:,3,zz) = foreground(:,:,zz)*128 + tip_dilate(:,:,zz)*255;
end
tifwrite(uint8(im), [file_name '\tip_tracking\' ind ]);



% % swc visiualization
% tt_write = 16;
% swc_path = [file_name 'swc\'  num2str(tt_write) '\'];
% mkdir(swc_path);
% tip_num = 0;
% for ii = 1:length(trajectories)
%     f_flag = 0;
%     for jj = 2:length(trajectories{ii})
%         if node_table(trajectories{ii}(jj),1) == tt_write && jj >= 8
%             f_flag = 1;
% %             break;
%         end
%     end
%     if f_flag == 1
%         tip_color = round(rand()*255 + 20);
%         tip_ind = trajectories{ii}(jj);
%         tip_num = tip_num + 1;
%         f = fopen([swc_path num2str(tip_num) '.swc'],'w');
%         node_ind = 1;
%         fprintf(f,'%d %d %d %d %d %d %d\n', node_ind, tip_color, node_table(tip_ind,3), 301- node_table(tip_ind,2), node_table(tip_ind,4), 3, -1);
%         fprintf(f,'%d %d %d %d %d %d %d\n', node_ind+1, tip_color, node_table(tip_ind,3), 301-node_table(tip_ind,2), node_table(tip_ind,4), 1, 1);
%         for kk = jj - 1 : -1 : 1
%           tip_ind = trajectories{ii}(kk);
%           node_ind = node_ind + 1;
%           fprintf(f,'%d %d %d %d %d %d %d\n', node_ind+1, tip_color, node_table(tip_ind,3), 301- node_table(tip_ind,2), node_table(tip_ind,4), 1, node_ind);
%         end
%         fclose(f);
% %         break
%     end   
% end


fore_tip = zeros(x,y,z);
track_tip = zeros(x,y,z);
tip1 = trajectories{ii_write}(jj_write);
fore_tip(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4)) = 1;
for jj =1:jj_write - 1
    
    tip1 = trajectories{ii_write}(jj);
    tt_list = node_table(trajectories{ii}(jj+1:end),1);
    track_tip(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4)) = 1;

end
ind = num2str(1000+tt_write);
ind = ind(2:4);
im = zeros(x, y, 3, z);
foreground = load([file_name 'foreground\' ind '.mat']);
foreground = foreground.foreground;
tip_dilate = imdilate(fore_tip,strel('sphere',3));
track_dilate = imdilate(track_tip,strel('sphere',1));
for zz = 1:z
    im(:,:,1,zz) = foreground(:,:,zz)*128;
    im(:,:,2,zz) = foreground(:,:,zz)*128 + tip_dilate(:,:,zz)*255 + track_dilate(:,:,zz)*255;
    im(:,:,3,zz) = foreground(:,:,zz)*128 + tip_dilate(:,:,zz)*255;
end
tifwrite(uint8(im), [file_name '\tip_tracking\' ind ]);

% fore_tip_r = zeros(x,y,z,t);
% fore_tip_g = zeros(x,y,z,t);
% fore_tip_b = zeros(x,y,z,t);
% track_tip_r = zeros(x,y,z,t);
% track_tip_g = zeros(x,y,z,t);
% track_tip_b = zeros(x,y,z,t);
% for ii = 1:length(trajectories)
% for jj =1:length(trajectories{ii}) - 1
%     
%     tip1 = trajectories{ii}(jj);
%     tt_list = node_table(trajectories{ii}(jj+1:end),1);
%     rgb = rand(3,1);
% %         line(:,:,:,:,tt) = draw3Dline(line(:,:,:,:,tt), node_table(tip1,2:4), node_table(tip2,2:4),3,[0 255 0]);
%     fore_tip_r(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),node_table(trajectories{ii}(jj),1)) = rgb(1);
%     fore_tip_g(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),node_table(trajectories{ii}(jj),1)) = rgb(2);
%     fore_tip_b(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),node_table(trajectories{ii}(jj),1)) = rgb(3);
%     track_tip_r(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),tt_list) = rgb(1);
%     track_tip_g(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),tt_list) = rgb(2);
%     track_tip_b(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),tt_list) = rgb(3);
% 
% end
% end
% 
% for tt = 1:t
%     ind = num2str(1000+tt);
%     ind = ind(2:4);
%     im = zeros(x, y, 3, z);
%     foreground = load([file_name 'foreground\' ind '.mat']);
%     foreground = foreground.foreground;
%     tip_dilate_r = imdilate(fore_tip_r(:,:,:,tt),strel('sphere',3));
%     tip_dilate_g = imdilate(fore_tip_g(:,:,:,tt),strel('sphere',3));
%     tip_dilate_b = imdilate(fore_tip_b(:,:,:,tt),strel('sphere',3));
%     track_dilate_r = imdilate(track_tip_r(:,:,:,tt),strel('sphere',1));
%     track_dilate_g = imdilate(track_tip_g(:,:,:,tt),strel('sphere',1));
%     track_dilate_b = imdilate(track_tip_b(:,:,:,tt),strel('sphere',1));
%     for zz = 1:z
%         im(:,:,1,zz) = foreground(:,:,zz)*128 + tip_dilate_r(:,:,zz)*255 + track_dilate_r(:,:,zz)*255;
%         im(:,:,2,zz) = foreground(:,:,zz)*128 + tip_dilate_g(:,:,zz)*255 + track_dilate_g(:,:,zz)*255;
%         im(:,:,3,zz) = foreground(:,:,zz)*128 + tip_dilate_b(:,:,zz)*255 + track_dilate_b(:,:,zz)*255;
%     end
%     tifwrite(uint8(im), [file_name '\tip_tracking\' ind ]);
% end