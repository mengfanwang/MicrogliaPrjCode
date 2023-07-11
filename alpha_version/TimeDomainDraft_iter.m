clc;clear;close all;

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
load([file_name '\data_c2']);
data = data2;
[x,y,z,t] = size(data);

%%% tracking
load([file_name '\tip_info.mat']);
addpath ./CINDA
v_max = 5;
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

phi_initial = zeros(x,y,z,3,t);
for ii = 1:length(trajectories)
    for jj = 1:length(trajectories{ii})-1
        sp = node_table(trajectories{ii}(jj),2:4);
        ep = node_table(trajectories{ii}(jj+1),2:4);
        tt = node_table(trajectories{ii}(jj),1);
%         phi_initial(sp(1),sp(2),sp(3),1,tt) = sp(1) - ep(1);
%         phi_initial(sp(1),sp(2),sp(3),2,tt) = sp(2) - ep(2);
%         phi_initial(sp(1),sp(2),sp(3),3,tt) = sp(3) - ep(3);
        min_x = max(1,sp(1)-2); max_x = min(x,sp(1)+2);
        min_y = max(1,sp(2)-2); max_y = min(y,sp(2)+2);
        min_z = max(1,sp(3)-1); max_z = min(z,sp(3)+1);
        phi_initial(min_x:max_x,min_y:max_y,min_z:max_z,1,tt) = sp(1) - ep(1);
        phi_initial(min_x:max_x,min_y:max_y,min_z:max_z,2,tt) = sp(2) - ep(2);
        phi_initial(min_x:max_x,min_y:max_y,min_z:max_z,3,tt) = sp(3) - ep(3);
    end
end


% % vornono digram initial
% arrow_list = cell(t-1,1);
% for ii = 1:length(trajectories)
%     for jj = 1:length(trajectories{ii})-1
%         sp = node_table(trajectories{ii}(jj),2:4);
%         ep = node_table(trajectories{ii}(jj+1),2:4);
%         tt = node_table(trajectories{ii}(jj),1);
%         arrow_list{tt} = [arrow_list{tt}; [sp ep]];
%     end
% end
% phi_initial = zeros(x,y,z,3,t);
% [x_ind,y_ind,z_ind] = ind2sub([x y z],1:x*y*z);
% for tt = 9:9
%     [~,min_ind] = min(dist(arrow_list{tt}(:,1:3),[x_ind; y_ind; z_ind]));
%     phi_initial_x = zeros(x,y,z);
%     phi_initial_y = zeros(x,y,z);
%     phi_initial_z = zeros(x,y,z);
%     for jj = 1:length(arrow_list{tt})
%         phi_initial_x(min_ind==jj) = arrow_list{tt}(jj,1) - arrow_list{tt}(jj,4);
%         phi_initial_y(min_ind==jj) = arrow_list{tt}(jj,2) - arrow_list{tt}(jj,5);
%         phi_initial_z(min_ind==jj) = arrow_list{tt}(jj,3) - arrow_list{tt}(jj,6);
%     end
%     phi_initial(:,:,:,1,tt) = phi_initial_x; 
%     phi_initial(:,:,:,2,tt) = phi_initial_y; 
%     phi_initial(:,:,:,3,tt) = phi_initial_z; 
% end

% tt = 9;
for tt = 9:9
[x,y,z,t] = size(data);
data1 = data(:,:,:,tt)/255;
data2 = data(:,:,:,tt+1)/255;

% [ux,uy,uz]=LK3D(data1,data2,3);
% v = sqrt(ux.^2 + uy.^2 + uz.^2);
% v(v>10) = 10;
% load([file_name '\result_gradientDescent\result_gradientDescent_' num2str(tt) '_0.01_0.2_1500.mat']);
% phi_original = zeros([size(ux) 3]);
% phi_original(:,:,:,1) = ux;
% phi_original(:,:,:,2) = uy;
% phi_original(:,:,:,3) = uz;
% load([file_name '\result_gradientDescent_iter1\result_gradientDescent_' num2str(tt) '_0.01_0.2_2000.mat']);
% ux = ux - phi_original(:,:,:,1);
% uy = uy - phi_original(:,:,:,2);
% uz = uz - phi_original(:,:,:,3);
ux = phi_initial(:,:,:,1,tt);
uy = phi_initial(:,:,:,2,tt); 
uz = phi_initial(:,:,:,3,tt);
%% error calculation
data1_pad = padarray(data1,[5 5 5],'replicate'); 
data2_pad = padarray(data2,[5 5 5],'replicate');
[x_mesh,y_mesh,z_mesh] = meshgrid(-4:x+5,-4:y+5,-4:z+5);
[x_q,y_q,z_q] = meshgrid(1:x,1:y,1:z);
data1_trans = interp3(x_mesh,y_mesh,z_mesh,data1_pad,x_q+uy,y_q+ux,z_q+uz);
% data1_trans = result;
fprintf('Original error:%f\n', sum((data1(:)-data2(:)).^2)/x/y/z);
fprintf('Similarity error:%f\n', sum((data1_trans(:)-data2(:)).^2)/x/y/z);
%
dux = zeros(size(ux));
duy = zeros(size(uy));
duz = zeros(size(uz));
% du/dx(x,y,z)
dux(1:end-1,:,:,1) = ux(2:end,:,:) - ux(1:end-1,:,:); 
dux(1:end-1,:,:,2) = uy(2:end,:,:) - uy(1:end-1,:,:);
dux(1:end-1,:,:,3) = uz(2:end,:,:) - uz(1:end-1,:,:);
% du/dy(x,y,z)
duy(:,1:end-1,:,1) = ux(:,2:end,:) - ux(:,1:end-1,:);  
duy(:,1:end-1,:,2) = uy(:,2:end,:) - uy(:,1:end-1,:);
duy(:,1:end-1,:,3) = uz(:,2:end,:) - uz(:,1:end-1,:);
% du/dz(x,y,z)
duz(:,:,1:end-1,1) = ux(:,:,2:end) - ux(:,:,1:end-1); 
duz(:,:,1:end-1,2) = uy(:,:,2:end) - uy(:,:,1:end-1);
duz(:,:,1:end-1,3) = uz(:,:,2:end) - uz(:,:,1:end-1);
    
dv = dux.^2 + duy.^2 + duz.^2;
num = (x-1)*y*z*3 + x*(y-1)*z*3 + x*y*(z-1)*3;
fprintf('Smoothness error:%f\n',sum(dv(:))/num);

% % save image
tifwrite(uint8(data1_trans*255),'temp_result');
tifwrite(uint8((data1-data2).^2*255),'temp_diff');
tifwrite(uint8((data1_trans-data2).^2*255),'temp_result_diff');
v = sqrt(ux.^2 + uy.^2 + uz.^2);
v = v.^2/max(v(:))*255;
% v(v<80) = 0;
tifwrite(uint8(v),'temp_field');
% load([file_name '\foreground_c2_reconnect\009' ]);
% v(foreground==0) = nan;
% writeColormap(v, 'temp_field');

% Enhance the quiver plot visually by downsizing vectors  
%   -f : downsizing factor
v = sqrt(ux.^2 + uy.^2 + uz.^2);
x = zeros(301,301); y = zeros(301,301);
[~, v_maxind] = max(v,[],3);
for ii = 1:301
    for jj = 1:301
        x(ii,jj) = ux(ii,jj,v_maxind(ii,jj));
        y(ii,jj) = uy(ii,jj,v_maxind(ii,jj));
    end
end
color = colormap;
color_quiver = ceil(max(v,[],3)*63/max(v(:)+1e-6)) + 1;

f=1;
x = x(1:f:301,1:f:301);
y = y(1:f:301,1:f:301);
[X,Y]=meshgrid(1:size(x,2),1:size(x,1));
scale = 4;
h = imshow(imresize(max(data1,[],3)*2/3,scale));hold on;
% h = imshow(imresize(max((data1-data2).^2*2/3,[],3),4));hold on;
for ii = 1:64
    ind = color_quiver == ii;
    quiver(X(ind)*scale,Y(ind)*scale,-y(ind)*scale,-x(ind)*scale,0,'Color',color(ii,:));
end
% quiver(X,Y,x,y); 
axis([1 size(x,2) 1 size(x,1)]);

ind = num2str(1000+tt);
ind = ind(2:4);
zoom(0.1);
% export_fig([file_name '\time_domain_diff\' ind '.png'], '-transparent', '-m4');
end

function writeColormap(distance_map, filename)
distance_map(distance_map>10) = 10;
distance_map = ceil(distance_map*63/max(distance_map(:)+1e-6)) + 1;
color = colormap;
[x,y,z] = size(distance_map);
map = zeros(x,y,3,z);
for ii = 1:x
    for jj = 1:y
        for kk = 1:z
            if ~isnan(distance_map(ii,jj,kk))
                map(ii,jj,1,kk) = color(distance_map(ii,jj,kk),1);
                map(ii,jj,2,kk) = color(distance_map(ii,jj,kk),2);
                map(ii,jj,3,kk) = color(distance_map(ii,jj,kk),3);
            end
        end
    end
end
tifwrite(uint8(map*255),filename);
end
    