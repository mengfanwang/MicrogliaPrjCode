clc;clear;close all;

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
load([file_name '\data_c2']);
data = data2;
[x,y,z,t] = size(data);

foreground = zeros(x,y,z,t);
for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    temp = load([file_name '\foreground_c2_iter1\' ind]);
    foreground(:,:,:,tt) = temp.foreground;
end

iter = 0;
[x_mesh,y_mesh,z_mesh] = meshgrid(-4:x+5,-4:y+5,-4:z+5);
[x_q,y_q,z_q] = meshgrid(1:x,1:y,1:z);

diff = 1;
while diff ~= 0
iter = iter + 1;
fore_new = zeros(x,y,z,t);
for tt = 2:t-1
    fore2 = foreground(:,:,:,tt);
    load([file_name '\result_gradientDescent\result_gradientDescent_' num2str(tt-1) '_0.01_0.2_1500.mat']);
    fore1 = foreground(:,:,:,tt-1);
    fore1_pad = padarray(fore1,[5 5 5],'replicate'); 
    fore1_forward = interp3(x_mesh,y_mesh,z_mesh,fore1_pad,x_q+uy,y_q+ux,z_q+uz);
    load([file_name '\result_gradientDescent_inverse\result_gradientDescent_' num2str(tt+1) '_0.01_0.2_1500.mat']);
    fore3 = foreground(:,:,:,tt+1);
    fore3_pad = padarray(fore3,[5 5 5],'replicate'); 
    fore3_backward = interp3(x_mesh,y_mesh,z_mesh,fore3_pad,x_q+uy,y_q+ux,z_q+uz);
    fore_temp = (fore1_forward + fore2 + fore3_backward)/3;

    ind = find(fore_temp>0.5/3 & fore_temp < 1/3);
    fore_temp(fore_temp<=0.5/3) = 0;
    fore_temp(fore_temp>=1/3) = 1;
    fore_temp(ind) = fore2(ind);
    fore_new(:,:,:,tt) = fore_temp;
end
diff = sum(abs(fore_new(:,:,:,2:t-1) - foreground(:,:,:,2:t-1)), 'all');
fprintf('Difference of iter %d: %d\n',iter,diff);
foreground(:,:,:,2:t-1) = fore_new(:,:,:,2:t-1);

end

%%
for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    tifwrite(uint8(foreground(:,:,:,tt)*255),[file_name '\foreground_c2_iter1_time_reconnect\' ind]);
end

