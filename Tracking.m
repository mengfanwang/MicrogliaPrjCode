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
v_max = 12;




% multiframe graph
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
trajectories(cellfun(@length,trajectories)<5) = [];

% swc visiualization
location = [140.911357 158.143121 139.799631 163.701754 160.088643 168.426593 161.756233 160.088643];
locx = mean(location(2:2:8));
locy = mean(location(1:2:7));
locz = 18.5;
tt_write = 13;
swc_path = [file_name 'swc\'  num2str(tt_write) '\'];
mkdir(swc_path);
tip_num = 0;
for ii = 1:length(trajectories)
    f_flag = 0;
    for jj = 2:length(trajectories{ii})
        if node_table(trajectories{ii}(jj),1) == tt_write && trajectories{ii}(jj) - trajectories{ii}(1) >= 5
            f_flag = 1;
%             break;
        end
    end
    tip1 = trajectories{ii}(jj);
    headYXZ = [node_table(tip1,2) node_table(tip1,3) node_table(tip1,4)];
    if f_flag == 1 && headYXZ(1) <= locx + 50 && headYXZ(1) >= locx - 50 && headYXZ(2) <= locy + 50 && headYXZ(2) >= locy - 50
        tip_color = round(rand()*255 + 20);
        tip_ind = trajectories{ii}(jj);
        tip_num = tip_num + 1;
        f = fopen([swc_path num2str(tip_num) '.swc'],'w');
        node_ind = 1;
        fprintf(f,'%d %d %d %d %d %d %d\n', node_ind, tip_color, node_table(tip_ind,3), 301- node_table(tip_ind,2), node_table(tip_ind,4), 3, -1);
        fprintf(f,'%d %d %d %d %d %d %d\n', node_ind+1, tip_color, node_table(tip_ind,3), 301-node_table(tip_ind,2), node_table(tip_ind,4), 1, 1);
        for kk = jj : -1 : 1
          tip_ind = trajectories{ii}(kk);
          node_ind = node_ind + 1;
          fprintf(f,'%d %d %d %d %d %d %d\n', node_ind+1, tip_color, node_table(tip_ind,3), 301- node_table(tip_ind,2), node_table(tip_ind,4), 1, node_ind);
        end
        fclose(f);
%         break
    end   
end

traj_maxnum = 0;
for ii = 1:t-1
    traj_maxnum = traj_maxnum + min(num_tip(ii),num_tip(ii+1));
end
disp('Recall:');
sum(cellfun(@length,trajectories)-1)/traj_maxnum

location = [140.911357 158.143121 139.799631 163.701754 160.088643 168.426593 161.756233 160.088643];
locx = mean(location(2:2:8));
locy = mean(location(1:2:7));
locz = 18.5;
forward = [];
for ii = 1:length(trajectories)
    for jj = 1:length(trajectories{ii})-1
        pix1 = trajectories{ii}(jj);
        pix2 = trajectories{ii}(jj+1);
        headYXZ = [node_table(pix1,2) node_table(pix1,3) node_table(pix1,4)];
        tailYXZ = [node_table(pix2,2) node_table(pix2,3) node_table(pix2,4)];
        if tailYXZ(1) <= locx + 50 && tailYXZ(1) >= locx - 50 && tailYXZ(2) <= locy + 50 && tailYXZ(2) >= locy - 50
            vec_a = tailYXZ - headYXZ;
            vec_b = [locx locy locz] - headYXZ;
            if norm(vec_a)~=0 && norm(vec_b)~=0
                forward = [forward vec_a*vec_b'/norm(vec_b)];
            end
        end
    end
end
mean(forward)


% draw 3D line
im = zeros(x, y, 3, z, t);
fore_tip_x = zeros(x,y,z,t);
fore_tip_y = zeros(x,y,z,t);
fore_tip_z = zeros(x,y,z,t);
for ii = 1:size(node_table,1)
    fore_tip_x(node_table(ii,2),node_table(ii,3),node_table(ii,4),node_table(ii,1)) = 255;
    fore_tip_y(node_table(ii,2),node_table(ii,3),node_table(ii,4),node_table(ii,1)) = 255;
    fore_tip_z(node_table(ii,2),node_table(ii,3),node_table(ii,4),node_table(ii,1)) = 255;
end
headYXZs = cell(t,1);
tailYXZs = cell(t,1);
colors = cell(t,1);
for ii = 1:length(trajectories)
    color = round(rand(1,3)*255);
    while sum(color) > 2.5*255 || sum(color) < 0.5*255
        color = round(rand(1,3)*255);
    end
    for jj = 1:length(trajectories{ii})-1
        pix1 = trajectories{ii}(jj);
        pix2 = trajectories{ii}(jj+1);
        headYXZ = [node_table(pix1,2) node_table(pix1,3) node_table(pix1,4)];
        tailYXZ = [node_table(pix2,2) node_table(pix2,3) node_table(pix2,4)];      
        
        fore_tip_x(node_table(pix1,2),node_table(pix1,3),node_table(pix1,4),node_table(pix1,1)) = color(1);
        fore_tip_y(node_table(pix1,2),node_table(pix1,3),node_table(pix1,4),node_table(pix1,1)) = color(2);
        fore_tip_z(node_table(pix1,2),node_table(pix1,3),node_table(pix1,4),node_table(pix1,1)) = color(3);
        fore_tip_x(node_table(pix2,2),node_table(pix2,3),node_table(pix2,4),node_table(pix2,1)) = color(1);
        fore_tip_y(node_table(pix2,2),node_table(pix2,3),node_table(pix2,4),node_table(pix2,1)) = color(2);
        fore_tip_z(node_table(pix2,2),node_table(pix2,3),node_table(pix2,4),node_table(pix2,1)) = color(3);
        
        for tt = node_table(pix2,1):node_table(trajectories{ii}(end),1)
%             im(:,:,:,:,tt) = draw3Dline(im(:,:,:,:,tt), headYXZ, tailYXZ, 1, [0 255 0]);
%             im_(:,:,:,:,tt) = draw3Dline(im_(:,:,:,:,tt), headYXZ*3-1, tailYXZ*3-1, 1, [0 255 0]);
              headYXZs{tt} = [headYXZs{tt}; headYXZ*3-1];
              tailYXZs{tt} = [tailYXZs{tt}; tailYXZ*3-1];
              colors{tt} = [colors{tt}; color];
        end
    end
end
for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    foreground = load([file_name 'foreground_c2\' ind '.mat']);
    foreground = foreground.foreground;
    tip_dilate_x = imdilate(fore_tip_x(:,:,:,tt),strel('sphere',2));
    tip_dilate_y = imdilate(fore_tip_y(:,:,:,tt),strel('sphere',2));
    tip_dilate_z = imdilate(fore_tip_z(:,:,:,tt),strel('sphere',2));
    for zz = 1:z
        im(:,:,1,zz,tt) = foreground(:,:,zz)*128.*(tip_dilate_x(:,:,zz)==0) + tip_dilate_x(:,:,zz);
        im(:,:,2,zz,tt) = foreground(:,:,zz)*128.*(tip_dilate_x(:,:,zz)==0) + tip_dilate_y(:,:,zz);
        im(:,:,3,zz,tt) = foreground(:,:,zz)*128.*(tip_dilate_x(:,:,zz)==0) + tip_dilate_z(:,:,zz);
    end
end
n = 3;
% im_ = im;
im_ = zeros(n*x, n*y, 3, n*z, t);
for tt = 1:t
    tt
    im_(:,:,:,:,tt) = imEnlarge(im(:,:,:,:,tt),n);
end
toc
for tt = 1:t
    tt
    ind = num2str(1000+tt);
    ind = ind(2:4);
    if tt > 1
        im_(:,:,:,:,tt) = draw3Dlines(im_(:,:,:,:,tt), headYXZs{tt}, tailYXZs{tt}, 1, colors{tt});
    end
%     tifwrite(uint8(im_(:,:,:,:,tt)), [file_name '\tip_tracking\' ind ]);
    tifwrite(uint8(im_(:,:,:,:,tt)), ['C:\Users\Mengfan Wang\Dropbox\temp_data\tip_tracking\' ind ]);
end
toc


% plot single_line
single_ind = 1694;
for ii = 1:length(trajectories)
    if ismember(1694, trajectories{ii})
        break
    end
end
fore_tip = zeros(x,y,z,t);
for ii = 1:length(trajectories)
for jj =1:length(trajectories{ii}) - 1
    
    tip1 = trajectories{ii}(jj);
    tt_list = node_table(trajectories{ii}(jj+1:end),1);
%         line(:,:,:,:,tt) = draw3Dline(line(:,:,:,:,tt), node_table(tip1,2:4), node_table(tip2,2:4),3,[0 255 0]);
    fore_tip(node_table(tip1,2),node_table(tip1,3),node_table(tip1,4),tt_list) = 1;

end
end
for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    im = tifread([file_name '\tip_connection\' ind '.tif']);
    tip_dilate = imdilate(fore_tip(:,:,:,tt),strel('sphere',1));
    for zz = 1:z
        im(:,:,2,zz) = im(:,:,2,zz) + tip_dilate(:,:,zz)*255;
    end
    tifwrite(uint8(im), [file_name '\tip_tracking\' ind ]);
end

function im_ = imEnlarge(im, n)
    % it can enlarge a 3d color image x times
    [x, y, c, z] = size(im);
    im_ = zeros(x*n, y*n, 3, z*n);
    for zz = 1:z
        im_(:,:,:,(zz-1)*n+1:zz*n) = repmat(imresize(im(:,:,:,zz),n),1,1,1,n);
    end
%     for xx = 1:x
%         for yy = 1:y
%             for zz = 1:z
%                 im_((xx-1)*n+1:xx*n, (yy-1)*n+1:yy*n, 1, (zz-1)*n+1:zz*n) = im(xx,yy,1,zz);
%                 im_((xx-1)*n+1:xx*n, (yy-1)*n+1:yy*n, 2, (zz-1)*n+1:zz*n) = im(xx,yy,2,zz);
%                 im_((xx-1)*n+1:xx*n, (yy-1)*n+1:yy*n, 3, (zz-1)*n+1:zz*n) = im(xx,yy,3,zz);
%             end
%         end
%     end
end