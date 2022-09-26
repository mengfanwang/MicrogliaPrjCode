clc;clear;close all;
dbstop if error
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';

%%% load gt %%%
traj_num = 75;
trajectories = cell(10,1);
t_ind = 0;
for tt = 1:traj_num
    ind = num2str(1000+tt);
    ind = ind(2:4);
    load([file_name 'gt\' ind '.mat']);
    if length(traj) >= 20
        t_ind = t_ind + 1;
        trajectories{t_ind} = traj;
    end
end


% load data
load([file_name '\data_c1.mat']);
load([file_name '\data_c2.mat']);
[x,y,z,t] = size(data2);
scale = 2;
size_new = round(scale*[1 1 1/0.461].*[x y z]);
x_n = size_new(1);
y_n = size_new(2);
z_n = size_new(3);
% data = imresize3(data,scale);
data_r = zeros(x_n, y_n, z_n, t);
data_g = zeros(x_n, y_n, z_n, t);
for tt = 1:t
    data_r(:,:,:,tt) = imresize3(data1(:,:,:,tt), size_new);
    data_g(:,:,:,tt) = imresize3(data2(:,:,:,tt), size_new);
end
data_c = zeros(x_n, y_n, 3, z_n, t);
for tt = 1:t
    for zz = 1:z_n
        data_c(:,:,1,zz,tt) = data_r(:,:,zz,tt);
        data_c(:,:,2,zz,tt) = data_g(:,:,zz,tt)/1.5;
% %         data_c(:,:,3,zz,tt) = imresize(data(:,:,zz,tt),scale)/2;
    end
end
% fprintf('Plot trajectories...');
% for ii = 1:length(trajectories)
%     fprintf('%d/%d\n',ii,length(trajectories));
%     color = (rand(1,3)/1.5 + 0.5/1.5)*255;
%     for jj = 2:length(trajectories{ii})   
%         sp = trajectories{ii}(jj-1,1:3);
%         ep = trajectories{ii}(jj,1:3);
%         tt = trajectories{ii}(jj,4);
% %         t_end = trajectories{ii}(end,4);
%         for ttt = tt:t
%             data_c(:,:,:,:,ttt) = draw3Dline(data_c(:,:,:,:,ttt), [sp(1)*scale sp(2)*scale sp(3)*scale], [ep(1)*scale ep(2)*scale ep(3)*scale], 1.5, color);
%         end
%     end
% end

% write4dTiffRGB(data_c, [file_name '\tip_tracking_gt\im.tiff']);
% data_z = max(data_c,[],4);
% data_z = reshape(data_z, [x*scale y*scale 3 t]);

%% load infomation from czi file
% basic information
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
czi_file = bfopen([file_name '.czi']);
ome = czi_file{1,4};
rx = double(ome.getPixelsPhysicalSizeX(0).value);
ry = double(ome.getPixelsPhysicalSizeY(0).value);
rz = double(ome.getPixelsPhysicalSizeZ(0).value);

x = czi_file{1,4}.getPixelsSizeX(0).getValue();
y = czi_file{1,4}.getPixelsSizeY(0).getValue();
z = czi_file{1,4}.getPixelsSizeZ(0).getValue();
c = czi_file{1,4}.getPixelsSizeC(0).getValue();


% ablation time detection
ome_data = char(czi_file{1,4}.dumpXML());
ome_data = regexp(ome_data,'(?<=DeltaT=")\d+.\d+(?=")','match');
time_minute = zeros(length(ome_data)/c/z,1);
for ii = 0:length(ome_data)/c/z-1
    time_minute(ii+1) = str2num(ome_data{ii*c*z+1})/60;
end
time_minute = time_minute - time_minute(1);

location = [140.911357 158.143121 139.799631 163.701754 160.088643 168.426593 161.756233 160.088643];

%%
im = zeros(x_n + z_n + 20, y_n + z_n + 20, 3, t);
% plot ablation site
% location = [158.143121 140.911357 163.701754 139.799631;...
%             163.701754 139.799631 168.426593 160.088643;...
%             168.426593 160.088643 160.088643 161.756233;...
%             160.088643 161.756233 158.143121 140.911357]*2;
location = [140.911357 158.143121 139.799631 163.701754;...
            139.799631 163.701754 160.088643 168.426593;...
            160.088643 168.426593 161.756233 160.088643;...
            161.756233 160.088643 140.911357 158.143121]*scale;
        
polygon = insertShape(zeros(size(im(:,:,:,tt))),'Line',location,'LineWidth',1);
for tt = 1:t
    tt
    % 3d data max projection
    for cc = 1:3
        data_x = max(data_c(:,:,cc,:,tt),[],1);
        data_x = reshape(data_x, [y_n z_n]);
        data_x = data_x';
        im( (x_n+21):end, 1:y_n, cc, tt) = data_x;
        
        data_y = max(data_c(:,:,cc,:,tt),[],2);
        data_y = reshape(data_y, [x_n z_n]);
        im(1:x_n, (y_n+21):end, cc, tt) = data_y;
        
        data_z = max(data_c(:,:,cc,:,tt),[],4);
        data_z = reshape(data_z, [x_n y_n]);
        im(1:x_n, 1:y_n, cc, tt) = data_z;
    end
    
    time_temp = num2str(100 + round(time_minute(tt)));
    time_temp = time_temp(2:3);
    time_add = insertText(zeros(size(im(:,:,:,tt))),[20 20], ['t = ' time_temp ' min'] ,'FontSize', 16, 'TextColor', 'w', 'BoxOpacity', 0);
    scale_add = insertShape(zeros(size(im(:,:,:,tt))),'Line',[510-43 590 510 590],'LineWidth',3, 'Color', 'w') + ...
        insertText(zeros(size(im(:,:,:,tt))),[458 560], ['10 um'] ,'FontSize', 16, 'TextColor', 'w', 'BoxOpacity', 0);
    im(:,:,:,tt) = im(:,:,:,tt) + (polygon + time_add + scale_add)*255;
end

% second way to plot trajs
im2 = zeros(size(im));
fprintf('Plot trajectories...');
for ii = 1:length(trajectories)
    fprintf('%d/%d\n',ii,length(trajectories));
    color = (rand(1,3)/1.5 + 0.5/1.5);
%     color = (rand(1,3));
    for jj = 2:length(trajectories{ii})   
        sp = trajectories{ii}(jj-1,1:3);
        ep = trajectories{ii}(jj,1:3);
        tt = trajectories{ii}(jj,4);
        for ttt = tt:t
%             data_c(:,:,:,:,ttt) = draw3Dline(data_c(:,:,:,:,ttt), [sp(1)*scale sp(2)*scale sp(3)*scale], [ep(1)*scale ep(2)*scale ep(3)*scale], 1.5, color);
            im2(:,:,:,ttt) = im2(:,:,:,ttt) + 255*insertShape(zeros(size(im(:,:,:,tt))),'Line',[sp(2) sp(1) ep(2) ep(1)].*scale,'LineWidth',2, 'Color', color);
            im2(:,:,:,ttt) = im2(:,:,:,ttt) + 255*insertShape(zeros(size(im(:,:,:,tt))),'Line',[sp(3)/0.461 sp(1) ep(3)/0.461 ep(1)].*scale + [y_n+20 0 y_n+20 0],'LineWidth',2, 'Color', color);
            im2(:,:,:,ttt) = im2(:,:,:,ttt) + 255*insertShape(zeros(size(im(:,:,:,tt))),'Line',[sp(2) sp(3)/0.461 ep(2) ep(3)/0.461].*scale + [0 x_n+20 0 x_n+20],'LineWidth',2, 'Color', color);
        end
    end
end
im(im2>0) = im2;
tifwrite(uint8(im), [file_name '\tip_tracking_gt\im3']);
% %% write video
% im(im>255) = 255;
% v = VideoWriter([file_name '\tip_tracking_gt\im.mp4'], 'MPEG-4');
% v.Quality = 95;
% v.FrameRate = 5;
% open(v);
% for tt = 1:t
%     writeVideo(v, im(:,:,:,tt)/255);
% end
% close(v);