clc;clear;close all;

data = load('D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\stabilization_c2\stabilization_data.mat');
data = data.data;
data = data/255;
[x,y,z,t] = size(data);
% img = zeros(x,y,3,z,t);
im_tip_r = zeros(size(data));
im_tip_g = zeros(size(data));
im_tip_b = zeros(size(data));
im_traj_r = zeros(size(data));
im_traj_g = zeros(size(data));
im_traj_b = zeros(size(data));

path = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\gt';
trajs = dir(path);
traj_img = zeros(size(data));
for ii = 3:length(trajs)
    traj = load([path '\' trajs(ii).name]);
    traj = traj.traj;
    color = rand(3,1)/2 + 0.5;
    t_end = traj(end,4);
    for jj = 1:size(traj,1)
        im_tip_r(traj(jj,1),traj(jj,2),traj(jj,3),traj(jj,4)) = color(1);
        im_tip_g(traj(jj,1),traj(jj,2),traj(jj,3),traj(jj,4)) = color(2);
        im_tip_b(traj(jj,1),traj(jj,2),traj(jj,3),traj(jj,4)) = color(3);
        im_traj_r(traj(jj,1),traj(jj,2),traj(jj,3),traj(jj,4):t_end) = color(1);
        im_traj_g(traj(jj,1),traj(jj,2),traj(jj,3),traj(jj,4):t_end) = color(2);
        im_traj_b(traj(jj,1),traj(jj,2),traj(jj,3),traj(jj,4):t_end) = color(3);
    end
end
for tt = 1:t
    dilate_tip_r = imdilate(im_tip_r(:,:,:,tt),strel('sphere',3));
    dilate_tip_g = imdilate(im_tip_g(:,:,:,tt),strel('sphere',3));
    dilate_tip_b = imdilate(im_tip_b(:,:,:,tt),strel('sphere',3));
    dilate_traj_r = imdilate(im_traj_r(:,:,:,tt),strel('sphere',1));
    dilate_traj_g = imdilate(im_traj_g(:,:,:,tt),strel('sphere',1));
    dilate_traj_b = imdilate(im_traj_b(:,:,:,tt),strel('sphere',1));
    im = zeros(x,y,3,z);
    for zz = 1:z
        im(:,:,1,zz) = data(:,:,zz,tt).*(dilate_tip_r(:,:,zz)==0).*(dilate_traj_r(:,:,zz)==0)/1.5 + dilate_tip_r(:,:,zz) + dilate_traj_r(:,:,zz);
        im(:,:,2,zz) = data(:,:,zz,tt).*(dilate_tip_r(:,:,zz)==0).*(dilate_traj_r(:,:,zz)==0)/1.5 + dilate_tip_g(:,:,zz) + dilate_traj_g(:,:,zz);
        im(:,:,3,zz) = data(:,:,zz,tt).*(dilate_tip_r(:,:,zz)==0).*(dilate_traj_r(:,:,zz)==0)/1.5 + dilate_tip_b(:,:,zz) + dilate_traj_b(:,:,zz);
    end
    ind = num2str(1000+tt);
    ind = ind(2:4);
    tifwrite(uint8(im*255),['D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\stabilization\' ind]);
end