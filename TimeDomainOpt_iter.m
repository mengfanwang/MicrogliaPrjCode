clc;clear;close all;
% hope to optimize the registration problem:
% min L = L_similarity + L_smooth

% test for multilayer interpolation
win_size = [2 2 2];
a = rand(1,3);
a = a + win_size'/2;
% [x,y] = meshgrid(0:1,0:1);
coef = 0;
for ii = 1:win_size(1)
    for jj = 1:win_size(2)
        for kk = 1:win_size(3)
            coef = coef + (win_size(1)/2-abs(a(1)-ii))*(win_size(2)/2-abs(a(2)-jj))*(win_size(3)/2-abs(a(3)-kk));
        end
    end
end
fprintf("Coefficient of interpolation:%f\n", coef);
%

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



% tt = 9;
for tt = 9:9
data1 = data(:,:,:,tt)/255;
data2 = data(:,:,:,tt+1)/255;

lambda = 0.01;
pad_size = [10 10 6];
data1_pad = padarray(data1,pad_size,'replicate'); 
data2_pad = padarray(data2,pad_size,'replicate');

% gt1 = zeros(x,y,z);
% gt2 = zeros(x,y,z);
% [x_ind,y_ind,z_ind] = ind2sub(size(data1),1:x*y*z);
% for ii = -win_size(1)/2+1:win_size(1)/2
%     for jj = -win_size(2)/2+1:win_size(2)/2
%         for kk = -win_size(3)/2+1:win_size(3)/2
%             coef_x = win_size(1)/2 - abs(ii);
%             coef_y = win_size(2)/2 - abs(jj);
%             coef_z = win_size(3)/2 - abs(kk);
%             ind_new = sub2ind(size(data1_pad),x_ind+ii+pad_size(1),y_ind+jj+pad_size(2),z_ind+kk+pad_size(3));
%             phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/coef;
%             gt1 = gt1 + reshape(phi_current_temp,[x y z]);
%             phi_current_temp = data2_pad(ind_new).*coef_x.*coef_y.*coef_z/coef;
%             gt2 = gt2 + reshape(phi_current_temp,[x y z]);
%         end
%     end
% end
gt2 = data2;

%                 ind_new = sub2ind(size(data1_pad),x_new+pad_size(1),y_new+pad_size(2),z_new+pad_size(3));
%                 phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/coef;
%                 data1_tran = data1_tran + reshape(phi_current_temp,[x y z]);
                
% phi_current = zeros(x,y,z,3);
phi_current = gpuArray(phi_initial(:,:,:,:,tt));
[x_ind,y_ind,z_ind] = ind2sub(size(data1),1:x*y*z);
iter = 0;
tic;
loss = gpuArray(zeros(100000,1));
time = gpuArray(zeros(100000,1));
step = 1e-6;
lr = 0.2;
decay = 0.9;
% while 1
%     iter = iter + 1;
for iter = 1:2000
    phi_previous = phi_current;
    phi_gradient = gpuArray(zeros(x,y,z,3));
    x_bias = reshape(phi_previous(:,:,:,1),[1 x*y*z]);
    y_bias = reshape(phi_previous(:,:,:,2),[1 x*y*z]);
    z_bias = reshape(phi_previous(:,:,:,3),[1 x*y*z]);
    
    x_new = x_ind + x_bias;
    y_new = y_ind + y_bias;
    z_new = z_ind + z_bias;
    data1_tran = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    
    x_new = x_new + step;
    data1_x_incre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    
    x_new = x_new - step;
    y_new = y_new + step;
    data1_y_incre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    
    y_new = y_new - step;
    z_new = z_new + step;
    data1_z_incre = interp3(data1_pad,y_new+pad_size(2),x_new+pad_size(1),z_new+pad_size(3));
    
    phi_gradient(:,:,:,1) = reshape((data1_tran-gt2(:)').*(data1_x_incre - data1_tran),[x y z])/step;
    phi_gradient(:,:,:,2) = reshape((data1_tran-gt2(:)').*(data1_y_incre - data1_tran),[x y z])/step;
    phi_gradient(:,:,:,3) = reshape((data1_tran-gt2(:)').*(data1_z_incre - data1_tran),[x y z])/step;
    data1_tran = reshape(data1_tran,[x y z]);
    
%     gradient_comparison(phi_gradient, phi_current, gt2, data1_pad);
%     smooth_comparison(phi_gradient, phi_previous);
    smooth_gradient = gpuArray(zeros(size(phi_previous)));
    smooth_gradient(1:end-1,:,:,:) = smooth_gradient(1:end-1,:,:,:) + phi_previous(1:end-1,:,:,:) - phi_previous(2:end,:,:,:);
    smooth_gradient(:,1:end-1,:,:) = smooth_gradient(:,1:end-1,:,:) + phi_previous(:,1:end-1,:,:) - phi_previous(:,2:end,:,:);
    smooth_gradient(:,:,1:end-1,:) = smooth_gradient(:,:,1:end-1,:) + phi_previous(:,:,1:end-1,:) - phi_previous(:,:,2:end,:);
    smooth_gradient(2:end,:,:,:) = smooth_gradient(2:end,:,:,:) + phi_previous(2:end,:,:,:) - phi_previous(1:end-1,:,:,:);
    smooth_gradient(:,2:end,:,:) = smooth_gradient(:,2:end,:,:) + phi_previous(:,2:end,:,:) - phi_previous(:,1:end-1,:,:);
    smooth_gradient(:,:,2:end,:) = smooth_gradient(:,:,2:end,:) + phi_previous(:,:,2:end,:) - phi_previous(:,:,1:end-1,:);
    phi_gradient = phi_gradient + lambda*smooth_gradient;
    
    phi_current = phi_current - lr*decay^floor(iter/200)*phi_gradient;
    phi_current(phi_current>5) = 5;
    phi_current(phi_current<-5) = -5;
%     new_ind = sub2ind(size(data1_pad),x_ind+x_bias+pad_size(1),y_ind+y_bias+pad_size(2),z_ind+z_bias+pad_size(3));
%     data1_tran = reshape(data1_pad(new_ind),[x y z]);
    mse = mean((data1_tran(:) - gt2(:)).^2);
    fprintf('Iteration %d\n Current error:%f Time:%f\n',iter, mse*1000000, toc);
    loss(iter) = mse;
    time(iter) = toc;
%     plot_quiver(phi_current,data1);hold on;
%     if toc > 3600
%         break;
%     end
end
loss = gather(loss(1:iter));
time = gather(time(1:iter));
phi_current = gather(phi_current);
ux = phi_current(:,:,:,1);
uy = phi_current(:,:,:,2);
uz = phi_current(:,:,:,3);
% save([file_name '\result_gradientDescent\result_gradientDescent_' num2str(win_size) '_' num2str(lambda) '_' num2str(lr) '_' num2str(iter)  '.mat'],'loss','time','ux','uy','uz');
save([file_name '\temp\5x5x3\result_gradientDescent_' num2str(tt) '_' num2str(lambda) '_' num2str(lr) '_' num2str(iter)  '.mat'],'loss','time','ux','uy','uz');
end
g = gpuDevice(1);
reset(g);   

function gradient_comparison(phi_gradient, phi_previous, gt2, data1_pad)
    [~, loc] = max(abs(phi_gradient(:)));
%     [loc_x,loc_y,loc_z,loc_d] = ind2sub(size(phi_gradient), loc);
    loc_x = 130;loc_y = 110;loc_z = 5; loc_d = 3;
    step = 1e-6;
    pad_size = [10 10 6];
    
    [x_ind,y_ind,z_ind] = ind2sub(size(gt2),1:numel(gt2));
    x_bias = reshape(phi_previous(:,:,:,1),[1 numel(gt2)]);
    y_bias = reshape(phi_previous(:,:,:,2),[1 numel(gt2)]);
    z_bias = reshape(phi_previous(:,:,:,3),[1 numel(gt2)]);
    tran_before = zeros(size(gt2));
    for ii = -1:2
        for jj = -1:2
            for kk = 0:1
                x_new = floor(x_ind + x_bias + ii);
                y_new = floor(y_ind + y_bias + jj);
                z_new = floor(z_ind + z_bias + kk);
                coef_x = 2 - abs(x_ind + x_bias - x_new);
                coef_y = 2 - abs(y_ind + y_bias - y_new);
                coef_z = 1 - abs(z_ind + z_bias - z_new);
                ind_new = sub2ind(size(data1_pad),x_new+pad_size(1),y_new+pad_size(2),z_new+pad_size(3));
                phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/16;
                tran_before = tran_before + reshape(phi_current_temp,size(gt2));
            end
        end
    end
    
    phi_previous(loc_x,loc_y,loc_z,loc_d) = phi_previous(loc_x,loc_y,loc_z,loc_d) + step;
    x_bias = reshape(phi_previous(:,:,:,1),[1 numel(gt2)]);
    y_bias = reshape(phi_previous(:,:,:,2),[1 numel(gt2)]);
    z_bias = reshape(phi_previous(:,:,:,3),[1 numel(gt2)]);
    tran_after = zeros(size(gt2));
    for ii = -1:2
        for jj = -1:2
            for kk = 0:1
                x_new = floor(x_ind + x_bias + ii);
                y_new = floor(y_ind + y_bias + jj);
                z_new = floor(z_ind + z_bias + kk);
                coef_x = 2 - abs(x_ind + x_bias - x_new);
                coef_y = 2 - abs(y_ind + y_bias - y_new);
                coef_z = 1 - abs(z_ind + z_bias - z_new);
                ind_new = sub2ind(size(data1_pad),x_new+pad_size(1),y_new+pad_size(2),z_new+pad_size(3));
                phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/16;
                tran_after = tran_after + reshape(phi_current_temp,size(gt2));
            end
        end
    end
    
    s1 = phi_gradient(loc_x,loc_y,loc_z,loc_d)*2/numel(gt2)
    s2 = mean((tran_after(:) - gt2(:)).^2 - (tran_before(:) - gt2(:)).^2)/step
%     s3 = 
    s1/s2
end
function smooth_comparison(phi_gradient, phi_previous)
    [~, loc] = max(abs(phi_gradient(:)));
%     [loc_x,loc_y,loc_z,loc_d] = ind2sub(size(phi_gradient), loc);
        loc_x = 1;loc_y = 110;loc_z = 1; loc_d = 1;
    step = 1e-6;
    
    % method 1
    ux1 = phi_previous(:,:,:,1);
    uy1 = phi_previous(:,:,:,2);
    uz1 = phi_previous(:,:,:,3);
    dux1 = zeros(size(ux1));
    duy1 = zeros(size(uy1));
    duz1 = zeros(size(uz1));
    dux1(1:end-1,:,:) = ux1(2:end,:,:) - ux1(1:end-1,:,:);
    duy1(:,1:end-1,:) = uy1(:,2:end,:) - uy1(:,1:end-1,:);
    duz1(:,:,1:end-1) = uz1(:,:,2:end) - uz1(:,:,1:end-1);
    dv1 = dux1.^2 + duy1.^2 + duz1.^2;
    
    phi_previous(loc_x,loc_y,loc_z,loc_d) = phi_previous(loc_x,loc_y,loc_z,loc_d) + step;
    ux2 = phi_previous(:,:,:,1);
    uy2 = phi_previous(:,:,:,2);
    uz2 = phi_previous(:,:,:,3);
    dux2 = zeros(size(ux2));
    duy2 = zeros(size(uy2));
    duz2 = zeros(size(uz2));
    dux2(1:end-1,:,:) = ux2(2:end,:,:) - ux2(1:end-1,:,:);
    duy2(:,1:end-1,:) = uy2(:,2:end,:) - uy2(:,1:end-1,:);
    duz2(:,:,1:end-1) = uz2(:,:,2:end) - uz2(:,:,1:end-1);
    dv2 = dux2.^2 + duy2.^2 + duz2.^2;
%     num = (x-1)*y*z + x*(y-1)*z + x*y*(z-1);
    s1 = (sum(dv2(:))-sum(dv1(:)))/step;
    
    % method 2
    dux3 = zeros([size(ux1) 3]);
    duy3 = zeros([size(uy1) 3]);
    duz3 = zeros([size(uz1) 3]);
    % du/dx(x,y,z)
    dux3(1:end-1,:,:,1) = ux1(2:end,:,:) - ux1(1:end-1,:,:); 
    dux3(1:end-1,:,:,2) = uy1(2:end,:,:) - uy1(1:end-1,:,:);
    dux3(1:end-1,:,:,3) = uz1(2:end,:,:) - uz1(1:end-1,:,:);
    % du/dy(x,y,z)
    duy3(:,1:end-1,:,1) = ux1(:,2:end,:) - ux1(:,1:end-1,:);  
    duy3(:,1:end-1,:,2) = uy1(:,2:end,:) - uy1(:,1:end-1,:);
    duy3(:,1:end-1,:,3) = uz1(:,2:end,:) - uz1(:,1:end-1,:);
    % du/dz(x,y,z)
    duz3(:,:,1:end-1,1) = ux1(:,:,2:end) - ux1(:,:,1:end-1); 
    duz3(:,:,1:end-1,2) = uy1(:,:,2:end) - uy1(:,:,1:end-1);
    duz3(:,:,1:end-1,3) = uz1(:,:,2:end) - uz1(:,:,1:end-1);
    dv3 = dux3.^2 + duy3.^2 + duz3.^2;
    
    dux4 = zeros([size(ux1) 3]);
    duy4 = zeros([size(uy1) 3]);
    duz4 = zeros([size(uz1) 3]);
    % du/dx(x,y,z)
    dux4(1:end-1,:,:,1) = ux2(2:end,:,:) - ux2(1:end-1,:,:); 
    dux4(1:end-1,:,:,2) = uy2(2:end,:,:) - uy2(1:end-1,:,:);
    dux4(1:end-1,:,:,3) = uz2(2:end,:,:) - uz2(1:end-1,:,:);
    % du/dy(x,y,z)
    duy4(:,1:end-1,:,1) = ux2(:,2:end,:) - ux2(:,1:end-1,:);  
    duy4(:,1:end-1,:,2) = uy2(:,2:end,:) - uy2(:,1:end-1,:);
    duy4(:,1:end-1,:,3) = uz2(:,2:end,:) - uz2(:,1:end-1,:);
    % du/dz(x,y,z)
    duz4(:,:,1:end-1,1) = ux2(:,:,2:end) - ux2(:,:,1:end-1); 
    duz4(:,:,1:end-1,2) = uy2(:,:,2:end) - uy2(:,:,1:end-1);
    duz4(:,:,1:end-1,3) = uz2(:,:,2:end) - uz2(:,:,1:end-1);
    dv4 = dux4.^2 + duy4.^2 + duz4.^2;
%     num = (x-1)*y*z*3 + x*(y-1)*z*3 + x*y*(z-1)*3
    s2 = (sum(dv4(:))-sum(dv3(:)))/step
%     s2 = sum(dux4(:).^2-dux3(:).^2)/step
%     s2/s1
%     s3 = 2*ux1(loc_x,loc_y,loc_z,loc_d) - ux1(loc_x+1,loc_y,loc_z,loc_d) - ux1(loc_x-1,loc_y,loc_z,loc_d) + ...
%          2*ux1(loc_x,loc_y,loc_z,loc_d) - ux1(loc_x,loc_y+1,loc_z,loc_d) - ux1(loc_x,loc_y-1,loc_z,loc_d) + ...
%          2*ux1(loc_x,loc_y,loc_z,loc_d) - ux1(loc_x,loc_y,loc_z+1,loc_d) - ux1(loc_x,loc_y,loc_z-1,loc_d);
     
  
    phi_previous(loc_x,loc_y,loc_z,loc_d) = phi_previous(loc_x,loc_y,loc_z,loc_d) - step;
    
    smooth_gradient = zeros(size(phi_previous));
    smooth_gradient(1:end-1,:,:,:) = smooth_gradient(1:end-1,:,:,:) + phi_previous(1:end-1,:,:,:) - phi_previous(2:end,:,:,:);
    smooth_gradient(:,1:end-1,:,:) = smooth_gradient(:,1:end-1,:,:) + phi_previous(:,1:end-1,:,:) - phi_previous(:,2:end,:,:);
    smooth_gradient(:,:,1:end-1,:) = smooth_gradient(:,:,1:end-1,:) + phi_previous(:,:,1:end-1,:) - phi_previous(:,:,2:end,:);
    smooth_gradient(2:end,:,:,:) = smooth_gradient(2:end,:,:,:) + phi_previous(2:end,:,:,:) - phi_previous(1:end-1,:,:,:);
    smooth_gradient(:,2:end,:,:) = smooth_gradient(:,2:end,:,:) + phi_previous(:,2:end,:,:) - phi_previous(:,1:end-1,:,:);
    smooth_gradient(:,:,2:end,:) = smooth_gradient(:,:,2:end,:) + phi_previous(:,:,2:end,:) - phi_previous(:,:,1:end-1,:);
    smooth_gradient(loc_x,loc_y,loc_z,loc_d)*2
end

function plot_quiver(phi_current, data1)     
    ux = phi_current(:,:,:,1);
    uy = phi_current(:,:,:,2);
    uz = phi_current(:,:,:,3);
    v = sqrt(ux.^2 + uy.^2 + uz.^2);
    max(v(:))
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
    f = 1;
    x = x(1:301,1:f:301);
    y = y(1:f:301,1:f:301);
    [X,Y]=meshgrid(1:size(x,2),1:size(x,1));
    imshow(max(data1,[],3)/2);hold on;
    for ii = 1:64
        ind = color_quiver == ii;
        quiver(X(ind),Y(ind),-y(ind),-x(ind),0,'Color',color(ii,:));
    end
    % quiver(X,Y,x,y); 
    axis([1 size(x,2) 1 size(x,1)]);
end

