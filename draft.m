clc;clear;close all;

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
load([file_name '\data_c2']);
% load([file_name '\stabilization_c2\stabilization_data.mat']);
data = data2;
[x,y,z,t] = size(data);

% for tt = 1:t
%     ind = num2str(1000+tt);
%     ind = ind(2:4);
%     load([file_name '\foreground_iter1\' ind]);
%     tifwrite(uint8(horzcat(data(:,:,:,tt).*foreground, data(:,:,:,tt).*(1-foreground))), [file_name '\foreground_iter1_segmentation\' ind]);
% end


% file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\temp';
% % data1 = load([file_name '\009_0.7_time.mat']);
% data1 = load('D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\foreground_c2\009.mat');
% data1 = data1.foreground;
% data2 = load([file_name '\009_0.7_time.mat']);
% % data2 = load('D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\foreground_c2\009.mat');
% components_cell = data2.components_cell;
% data2 = data2.foreground;

% [x,y,z] = size(data1);
for tt = 1:1
foreground = zeros(x,y,3,z);
ind = num2str(1000+tt);
ind = ind(2:4);
% data1 = load([file_name '\foreground_c2_iter1\' ind]);
% data1 = data1.foreground;
data1 = logical(tifread([file_name '\foreground_iter1_simulation\001_ori.tif']));
% data2 = load([file_name '\foreground_c2_reconnect_iter1\' ind]);
% data2 = data2.foreground;
data2 = logical(tifread([file_name '\foreground_iter1_simulation\001.tif']));
for zz = 1:z
    foreground(:,:,1,zz) = data2(:,:,zz); 
    foreground(:,:,2,zz) = data1(:,:,zz) & data2(:,:,zz);
    foreground(:,:,3,zz) = data1(:,:,zz) & data2(:,:,zz);
end
tifwrite(uint8(foreground*255),[file_name '\foreground_iter1_simulation\comp']);
end