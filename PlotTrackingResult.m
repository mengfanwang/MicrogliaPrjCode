clc;clear;close all;
dbstop if error
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\redetection_info_v2';

%%% tracking
load(fullfile(file_name, 'tipDetection', 'tip_info.mat'));
load(['D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\data_c2.mat']);
data = data2(:,:,:,2:end-1);
% data = data2;
% read data
% data = zeros(301,301,41,33);
% for tt = 1:33
%     ind2 = num2str(1001+tt);
%     ind2 = ind2(2:4);
%     data_temp = tifread(fullfile('D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\time_domain_data_after_VST', [ind2 '.tif']));
%     data(:,:,:,tt) = double(data_temp);
% end
[x,y,z,t] = size(data);
addpath ../CINDA
v_max = 10;
num_tip = cellfun(@length,xCoord);
% cum_tip = cumsum(num_tip);
detection_arcs = zeros(sum(num_tip),4);
node_table = zeros(sum(num_tip),4);
node_matrix = zeros(size(data));
transition_arcs = zeros(sum(num_tip(1:t-1).*num_tip(2:t)),3);
transition_ind = 0;
for tt = 1:t
    for jj = 1:num_tip(tt)
        ind = sum(num_tip(1:tt-1))+jj;
        detection_arcs(ind,:) = [ind v_max/2 v_max/2 -v_max];
        node_table(ind,:) = [tt xCoord{tt}(jj) yCoord{tt}(jj) zCoord{tt}(jj)];
        node_matrix(xCoord{tt}(jj),yCoord{tt}(jj),zCoord{tt}(jj),tt) = ind;
        if tt ~= t
            for kk = 1:num_tip(tt+1)
                transition_ind = transition_ind + 1;
                transition_arcs(transition_ind,:) = [ind sum(num_tip(1:tt))+kk ...
                    sqrt((xCoord{tt}(jj)-xCoord{tt+1}(kk))^2 + (yCoord{tt}(jj)-yCoord{tt+1}(kk))^2 + (zCoord{tt}(jj)-zCoord{tt+1}(kk))^2)];
            end
        end
    end
end

[trajectories, costs] = mcc4mot(detection_arcs,transition_arcs);
trajectories(cellfun(@length,trajectories)<5) = [];

%% plot in 3d
scale = 1;
data_c = zeros(x*scale,y*scale,3,z,t);
for tt = 1:t
    for zz = 1:z
        data_c(:,:,1,zz,tt) = 0;
        data_c(:,:,2,zz,tt) = imresize(data(:,:,zz,tt),scale);
        data_c(:,:,3,zz,tt) = 0;
    end
end
data_c = data_c * 0.7;
fprintf('Plot trajectories...');
for ii = 1:length(trajectories)
    fprintf('%d/%d\n',ii,length(trajectories));
    color = (rand(1,3)/1.5+0.5/1.5)*255;
    for jj = 2:length(trajectories{ii})
        sp = node_table(trajectories{ii}(jj-1),2:4);
        ep = node_table(trajectories{ii}(jj),2:4);
        tt = node_table(trajectories{ii}(jj),1);
        for ttt = tt:t
            data_c(:,:,:,:,ttt) = draw3Dline(data_c(:,:,:,:,ttt), [sp(1)*scale sp(2)*scale sp(3)], [ep(1)*scale ep(2)*scale ep(3)], 1.5, color);
        end
    end
end
mkdir(fullfile(file_name, 'tracking'));
% write4dTiffRGB(data_c, [file_name '\tip_tracking_gt\im.tiff']);
for tt = 1:t
%     ind = num2str(1000+tt);
%     ind = ind(2:4);
    tifwrite(uint8(data_c(:,:,:,:,tt)),fullfile(file_name, 'tracking', num2str(tt)));
end

%% plot in 2d
scale = 2;
data_c = zeros(x*scale,y*scale,3,t);
for tt = 1:t
    data_c(:,:,1,tt) = 0;
    data_c(:,:,2,tt) = imresize(max(data(:,:,:,tt),[],3),scale);
    data_c(:,:,3,tt) = 0;
end
data_c = data_c/255;
for ii = 1:length(trajectories)
    fprintf('%d/%d\n',ii,length(trajectories));
    color = (rand(1,3)/1.5+0.5/1.5);
    for jj = 2:length(trajectories{ii})
        sp = node_table(trajectories{ii}(jj-1),2:4);
        ep = node_table(trajectories{ii}(jj),2:4);
        tt = node_table(trajectories{ii}(jj),1);
        for ttt = tt:node_table(trajectories{ii}(end),1)
%             data_c(:,:,:,:,ttt) = draw3Dline(data_c(:,:,:,:,ttt), [sp(1)*scale sp(2)*scale sp(3)], [ep(1)*scale ep(2)*scale ep(3)], 1.5, color);
            data_c(:,:,:,ttt) = insertShape(data_c(:,:,:,ttt),'Line',[sp(2) sp(1) ep(2) ep(1)]*scale,'LineWidth',3, 'Color', color);
        end
    end
end
mkdir(fullfile(file_name, 'tracking'));
tifwrite(uint8(data_c*255),fullfile(file_name, 'tracking', '2d_result'));

%% ploting tracked balls
loc_ball = zeros(x,y,z,t);
for ii = 1:length(trajectories)
    fprintf('%d/%d\n',ii,length(trajectories));
    for jj = 1:length(trajectories{ii})
        ep = node_table(trajectories{ii}(jj),2:4);
        tt = node_table(trajectories{ii}(jj),1);
        loc_ball(ep(1),ep(2),ep(3),tt) = ii;
    end
end
1
for tt = 1:t
    loc_ball(:,:,:,tt) = imdilate(loc_ball(:,:,:,tt),strel('sphere',3));
end
2
color = (rand(length(trajectories),3)/2+0.5);
im_ball = zeros(x,y,3,z,t);
for tt = 1:t
    for zz = 1:z
        im_ball(:,:,1,zz,tt) = data(:,:,zz,tt)*0.6;
        im_ball(:,:,2,zz,tt) = data(:,:,zz,tt)*0.6;
        im_ball(:,:,3,zz,tt) = data(:,:,zz,tt)*0.6;
    end
end
3
for ii = 1:length(trajectories)
    [x_ind, y_ind, z_ind, t_ind] = ind2sub(size(loc_ball), find(loc_ball == ii));
    for jj = 1:length(x_ind)
        im_ball(x_ind(jj), y_ind(jj), 1, z_ind(jj), t_ind(jj)) = max(im_ball(x_ind(jj), y_ind(jj), 1, z_ind(jj), t_ind(jj)), color(ii,1)*255);
        im_ball(x_ind(jj), y_ind(jj), 2, z_ind(jj), t_ind(jj)) = max(im_ball(x_ind(jj), y_ind(jj), 2, z_ind(jj), t_ind(jj)), color(ii,2)*255);
        im_ball(x_ind(jj), y_ind(jj), 3, z_ind(jj), t_ind(jj)) = max(im_ball(x_ind(jj), y_ind(jj), 3, z_ind(jj), t_ind(jj)), color(ii,3)*255);
    end
end
for tt = 1:t
    tifwrite(uint8(im_ball(:,:,:,:,tt)),fullfile(file_name, 'ball_5', num2str(tt)));
end
