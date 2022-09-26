clc;clear;close all;
dbstop if error
% hope to optimize the registration problem:
% min L = L_similarity + L_smooth
win_size = [2 2 2];coef = 1;
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
load([file_name '\data_c2']);
data = data2;
[x,y,z,t] = size(data);

tt = 9;
data1 = data(:,:,:,tt)/255;
data2 = data(:,:,:,tt+1)/255;

lambda = 0.02;
pad_size = [10 10 6];
data1_pad = padarray(data1,pad_size,'replicate'); 
data2_pad = padarray(data2,pad_size,'replicate');

gt1 = data1; gt2 = data2;

phi_current = zeros(x,y,z,3);
[x_ind,y_ind,z_ind] = ind2sub(size(data1),1:x*y*z);
iter = 0;
tic;
loss = zeros(100000,1);
time = zeros(100000,1);
step = 1e-6;
lr = 0.2;
decay = 0.9;
while 1
    iter = iter + 1;
    phi_previous = phi_current;

    [phi_gradient, data1_tran] = GradientDescent1(phi_previous,data1_pad,pad_size,gt2);
    [phi_gradient2, data1_tran2] = GradientDescent2(phi_previous,data1_pad,pad_size,gt2);
    phi_gradient2(phi_gradient2>5) = 5;
    phi_gradient2(phi_gradient2<-5) = -5;

    
    phi_current = phi_current - lr*decay^floor(iter/200)*phi_gradient;
    phi_current(phi_current>5) = 5;
    phi_current(phi_current<-5) = -5;



    
    
    mse1 = mean((phi_gradient(:) - phi_gradient2(:)).^2)/ mean(phi_gradient(:).^2);
    fprintf('Iteration %d\n Current error:%f Time:%f\n',iter, mse1, toc);
%     loss(iter) = mse;
%     time(iter) = toc;
end

function [phi_gradient, data1_tran] = GradientDescent1(phi_previous,data1_pad,pad_size,gt2)
    [x,y,z] = size(gt2);
    [x_ind,y_ind,z_ind] = ind2sub(size(gt2),1:x*y*z);
    phi_gradient = zeros(x,y,z,3);
    step = 1e-6;
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
end

function [phi_gradient, data1_tran] = GradientDescent2(phi_previous,data1_pad,pad_size,gt2)
    [x,y,z] = size(gt2);
    [x_ind,y_ind,z_ind] = ind2sub(size(gt2),1:x*y*z);
    phi_gradient = zeros(x,y,z,3);
    step = 1e-6;
    win_size = [2 2 2];
    coef = 1;
    x_bias = reshape(phi_previous(:,:,:,1),[1 x*y*z]);
    y_bias = reshape(phi_previous(:,:,:,2),[1 x*y*z]);
    z_bias = reshape(phi_previous(:,:,:,3),[1 x*y*z]);

    data1_tran = zeros(1,numel(gt2));
    data1_x_incre = zeros(1,numel(gt2));
    data1_y_incre = zeros(1,numel(gt2));
    data1_z_incre = zeros(1,numel(gt2));
    for ii = -win_size(1)/2+1:win_size(1)/2
        for jj = -win_size(2)/2+1:win_size(2)/2
            for kk = -win_size(3)/2+1:win_size(3)/2
                x_new = floor(x_ind + x_bias + ii);
                y_new = floor(y_ind + y_bias + jj);
                z_new = floor(z_ind + z_bias + kk);
                coef_x = win_size(1)/2 - abs(x_ind + x_bias - x_new);
                coef_y = win_size(2)/2 - abs(y_ind + y_bias - y_new);
                coef_z = win_size(3)/2 - abs(z_ind + z_bias - z_new);
                ind_new = sub2ind(size(data1_pad),x_new+pad_size(1),y_new+pad_size(2),z_new+pad_size(3));
                phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/coef;
                data1_tran = data1_tran + phi_current_temp;
                
                x_bias = x_bias + step;
                x_new = floor(x_ind + x_bias + ii);
%                 y_new = floor(y_ind + y_bias + jj);
%                 z_new = floor(z_ind + z_bias + kk);
                coef_x = win_size(1)/2 - abs(x_ind + x_bias - x_new);
%                 coef_y = win_size(2)/2 - abs(y_ind + y_bias - y_new);
%                 coef_z = win_size(3)/2 - abs(z_ind + z_bias - z_new);
                ind_new = sub2ind(size(data1_pad),x_new+pad_size(1),y_new+pad_size(2),z_new+pad_size(3));
                phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/coef;
                data1_x_incre = data1_x_incre + phi_current_temp;
                
                x_bias = x_bias - step;
                y_bias = y_bias + step;
                x_new = floor(x_ind + x_bias + ii);
                y_new = floor(y_ind + y_bias + jj);
%                 z_new = floor(z_ind + z_bias + kk);
                coef_x = win_size(1)/2 - abs(x_ind + x_bias - x_new);
                coef_y = win_size(2)/2 - abs(y_ind + y_bias - y_new);
%                 coef_z = win_size(3)/2 - abs(z_ind + z_bias - z_new);
                ind_new = sub2ind(size(data1_pad),x_new+pad_size(1),y_new+pad_size(2),z_new+pad_size(3));
                phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/coef;
                data1_y_incre = data1_y_incre + phi_current_temp;
                
                y_bias = y_bias - step;
                z_bias = z_bias + step;
%                 x_new = floor(x_ind + x_bias + ii);
                y_new = floor(y_ind + y_bias + jj);
                z_new = floor(z_ind + z_bias + kk);
%                 coef_x = win_size(1)/2 - abs(x_ind + x_bias - x_new);
                coef_y = win_size(2)/2 - abs(y_ind + y_bias - y_new);
                coef_z = win_size(3)/2 - abs(z_ind + z_bias - z_new);
                ind_new = sub2ind(size(data1_pad),x_new+pad_size(1),y_new+pad_size(2),z_new+pad_size(3));
                phi_current_temp = data1_pad(ind_new).*coef_x.*coef_y.*coef_z/coef;
                data1_z_incre = data1_z_incre + phi_current_temp;
                z_bias = z_bias - step;
            end
        end
    end
    
    phi_gradient(:,:,:,1) = reshape((data1_tran-gt2(:)').*(data1_x_incre - data1_tran),[x y z])/step;
    phi_gradient(:,:,:,2) = reshape((data1_tran-gt2(:)').*(data1_y_incre - data1_tran),[x y z])/step;
    phi_gradient(:,:,:,3) = reshape((data1_tran-gt2(:)').*(data1_z_incre - data1_tran),[x y z])/step;

    data1_tran = reshape(data1_tran,[x y z]);
    
end