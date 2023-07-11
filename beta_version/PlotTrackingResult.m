clc;clear;close all;
dbstop if error
file_name = 'D:\dropbox\Modify Series Data-3nd\Priority 3 files\HA-100920-slice2-cortex-vessel-Modify Series\';

%%% tracking
load([file_name 'registration_c2.mat']);
data_all = data2;
data = data2;
[x,y,z,t] = size(data_all);

%% tracking
load([file_name 'tip_info.mat']);
addpath ./CINDA
v_max = 10;
num_tip = cellfun(@length,xCoord);
% cum_tip = cumsum(num_tip);
detection_arcs = zeros(sum(num_tip),4);
node_table = zeros(sum(num_tip),4);
node_matrix = zeros(size(data_all));
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
trajectories(cellfun(@length,trajectories)<4) = [];

scale = 1;
data_c = zeros(x*scale,y*scale,3,z,t);
for tt = 1:t
    for zz = 1:z
        data_c(:,:,1,zz,tt) = imresize(data(:,:,zz,tt),scale)/2;
        data_c(:,:,2,zz,tt) = imresize(data(:,:,zz,tt),scale)/2;
        data_c(:,:,3,zz,tt) = imresize(data(:,:,zz,tt),scale)/2;
    end
end
fprintf('Plot trajectories...');
for ii = 31:31   %1:length(trajectories)
    fprintf('%d/%d\n',ii,length(trajectories));
    color = (rand(1,3)/1.5+0.5/1.5)*255;
    for jj = 2:length(trajectories{ii})
        sp = node_table(trajectories{ii}(jj-1),2:4);
        ep = node_table(trajectories{ii}(jj),2:4);
        tt = node_table(trajectories{ii}(jj),1);
%         data_c(:,:,:,:,tt) = draw3Dline(data_c(:,:,:,:,tt), [sp(1)*scale sp(2)*scale sp(3)], [ep(1)*scale ep(2)*scale ep(3)], 1.5, color);
        for ttt = tt:node_table(trajectories{ii}(end),1) %t
            data_c(:,:,:,:,ttt) = draw3Dline(data_c(:,:,:,:,ttt), [sp(1)*scale sp(2)*scale sp(3)], [ep(1)*scale ep(2)*scale ep(3)], 1.5, color);
        end
    end
end
mkdir([file_name '\tracking']);
for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    tifwrite(uint8(data_c(:,:,:,:,tt)),[file_name '\tracking\' num2str(ind)]);
end