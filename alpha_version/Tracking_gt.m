clc;clear;close all;
dbstop if error

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';
load([file_name 'stabilization_c2\stabilization_data.mat']);
data_all = data;
[x, y, z, t] = size(data_all);
addpath ./CINDA

trajs_name = dir([file_name '\gt']);
trajs = cell(length(trajs_name)-2,1);
for ii = 3:length(trajs_name)
    traj = load([file_name '\gt\' trajs_name(ii).name]);
    trajs{ii-2} = traj.traj;
end
tip_ind = 0;
coord = cell(35,1);
for ii = 1:length(trajs)
    for jj = 1:size(trajs{ii},1)
        tip_ind = tip_ind + 1;
        trajs{ii}(jj,5) = tip_ind;
        coord{trajs{ii}(jj,4)} = [coord{trajs{ii}(jj,4)}; trajs{ii}(jj,1:3) trajs{ii}(jj,5)];
    end
end

%%% parameter setting %%%
v_max = 10;




% multiframe graph
num_tip = cellfun(@numel,coord)/4;
detection_arcs = zeros(sum(num_tip),4);
node_table = zeros(sum(num_tip),4);
node_matrix = zeros(size(data_all));
transition_arcs = zeros(sum(num_tip(1:t-1).*num_tip(2:t)),3);
transition_ind = 0;
for tt = 1:t
    for jj = 1:num_tip(tt)
        ind = coord{tt}(jj,4);
        detection_arcs(ind,:) = [ind v_max/2 v_max/2 -v_max];
        node_table(ind,:) = [tt coord{tt}(jj,1) coord{tt}(jj,2) coord{tt}(jj,3)];
        node_matrix(coord{tt}(jj,1),coord{tt}(jj,2),coord{tt}(jj,3),tt) = ind;
        if tt ~= t
            for kk = 1:num_tip(tt+1)
                transition_ind = transition_ind + 1;
                transition_arcs(transition_ind,:) = [ind coord{tt+1}(kk,4) ...
                    sqrt((coord{tt}(jj,1)-coord{tt+1}(kk,1))^2 + (coord{tt}(jj,2)-coord{tt+1}(kk,2))^2 + (coord{tt}(jj,3)-coord{tt+1}(kk,3))^2)];
            end
        end
    end
end

[trajectories, costs] = mcc4mot(detection_arcs,transition_arcs);
trajectories(cellfun(@length,trajectories)<=1) = [];

trajectories_gt = zeros(sum(cellfun(@numel,trajs)/5-1),2);
traj_ind = 0;
for ii = 1:length(trajs)
    for jj = 1:size(trajs{ii},1)-1
        traj_ind = traj_ind + 1;
        trajectories_gt(traj_ind,:) = [trajs{ii}(jj,5) trajs{ii}(jj+1,5)];
    end
end

trajectories_detect = zeros(sum(cellfun(@length,trajectories)-1),2);
traj_ind = 0;
for ii = 1:length(trajectories)
    for jj = 1:size(trajectories{ii},1)-1
        traj_ind = traj_ind + 1;
        trajectories_detect(traj_ind,:) = [trajectories{ii}(jj) trajectories{ii}(jj+1)];
    end
end

correct = 0;
error = [];
for ii = 1:size(trajectories_detect,1)
    if ismember(trajectories_detect(ii,:),trajectories_gt,'row')
        correct = correct + 1;
    else
        error = [error; trajectories_detect(ii,:)];
    end
end
correct

for ii = 1:length(trajectories)
    ii
    flag = 0;
    for jj = 1:length(trajs)
        if length(trajectories{ii}) == length(trajs{jj}(:,5))
        if all(trajectories{ii} == trajs{jj}(:,5))
            flag = 1;
        end
        end
    end
    if flag == 0
        a = 1;
    end
end

% write wrong connection
% tip = zeros(x,y,z,t);
% for ii = 1:size(error,1)
%     for jj = 1:2
%         node = node_table(error(ii,jj),:);
%         tip(node(2),node(3),node(4),node(1)) = 1;
%     end
% end
% for tt = 1:t
%     ind = num2str(1000+tt);
%     ind = ind(2:4);
%     im = tip(:,:,:,tt);
%     im = imdilate(im,strel('sphere',1));
%     im = im + (1-im).*data_all(:,:,:,tt)/255/2;
%     tifwrite(uint8(im*255),[file_name 'error_check\' ind]);
% end





