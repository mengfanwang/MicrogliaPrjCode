clc;clear;close all;

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
load([file_name '\data_c2']);
data = data2;
[x,y,z,t] = size(data);

% tt = 9;
for tt = 34:34
[x,y,z,t] = size(data);
data1 = data(:,:,:,tt)/255;
data2 = data(:,:,:,tt+1)/255;

% [ux,uy,uz]=LK3D(data1,data2,3);
% v = sqrt(ux.^2 + uy.^2 + uz.^2);
% v(v>10) = 10;'
% load([file_name '\temp\VoxelMorph_9\result_voxelMorph_20_1000_0.001_NCC.mat']);
load([file_name '\result_gradientDescent_inverse\result_gradientDescent_34_0.01_0.2_1500.mat']);
% load([file_name '\temp\5x5x5\result_gradientDescent_' num2str(tt) '_0.01_0.2_2000.mat']);

% load([file_name '\drift_c1.mat']);
% drift = -shifts1(tt).shifts + shifts1(tt+1).shifts;
% ux = ones(size(data1))*drift(1);
% uy = ones(size(data1))*drift(2);
% uz = ones(size(data1))*drift(3);
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
% tifwrite(uint8(v),'temp_field');
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
color_quiver = ceil(sqrt(max(v,[],3))*63/sqrt(max(v(:))+1e-6)) + 1; % sqrt: change vector color

f=1;
x = x(1:f:301,1:f:301);
y = y(1:f:301,1:f:301);
[X,Y]=meshgrid(1:size(x,2),1:size(x,1));
scale = 2;
% h = imshow(imresize(max(data1,[],3)/2,scale));hold on;
% h = imshow(imresize(max((data1-data2).^2*2/3,[],3),4));hold on;
h = imshow(zeros(size(imresize(max(data1,[],3)/2,scale)))); hold on;
for ii = 1:64
    ind = color_quiver == ii;
    quiver(X(ind)*scale,Y(ind)*scale,-y(ind)*scale*2,-x(ind)*scale*2,0,'Color',color(ii,:));
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

