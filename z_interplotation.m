clc;clear;
file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
load([file_name '\data_c2']);
data = data2;
data = double(data);

res = [0.4613 0.4613 1];
[x_size, y_size, z_size, t_size] = size(data);

% data_new = zeros(x_size,y_size,2*z_size-1,t_size);
% for tt = 1:1
%     for zz = 1:z_size-1
%         data_new(:,:,2*zz-1,tt) = data(:,:,zz,tt);
%         data_new(:,:,2*zz  ,tt) = (data(:,:,zz,tt) + data(:,:,zz+1,tt))/2;
%     end
%     data_new(:,:,end,tt) = data(:,:,end,tt);
% end

[x_mesh, y_mesh, z_mesh] = meshgrid(0:res(1):(x_size-1)*res(1),...
                                    0:res(2):(y_size-1)*res(2), 0:res(3):(z_size-1)*res(3));
[x_new,  y_new,  z_new ] = meshgrid(0:res(1):(x_size-1)*res(1),...
                                    0:res(2):(y_size-1)*res(2), 0:res(1):(z_size-1)*res(3));
for tt = 1:1
    data_new(:,:,:,tt) = interp3(x_mesh, y_mesh, z_mesh, data(:,:,:,tt), x_new, y_new, z_new, 'spline');
end



ind = [file_name '\data_c2\z_' ];
tifwrite(uint8(data_new(:,:,:,1)),ind);
    