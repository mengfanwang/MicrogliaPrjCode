clc;clear;close all;
dbstop if error

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';

% %% plot ablation site
% location = [140.911357 158.143121 139.799631 163.701754 160.088643 168.426593 161.756233 160.088643];
% load([file_name '\data_c1.mat']);
% [x, y, z, t] = size(data1);
% im = zeros(x,y,t);
% for tt = 1:t
%     ind = num2str(1000+tt);
%     ind = ind(2:4);
%     data = data1(:,:,:,tt);
% 
%     data = max(data, [], 3);
% %     data = data/1.5;
%     for ii = 1:4
%         locx = round(location(ii*2));
%         locy = round(location(ii*2-1));
%         data(locx-1:locx+1,locy-1:locy+1) = 255;
%     end
%    	im(:,:,tt) = data;
% end
% tifwrite(uint8(im),[file_name 'temp']);

% plot abalation site in 3D
location = [140.911357 158.143121 139.799631 163.701754 160.088643 168.426593 161.756233 160.088643];
load([file_name '\foreground_c1.mat']);
[x, y, z] = size(foreground);
a(1) = round(location(2));
a(2) = round(location(1));
b(1) = round(location(4));
b(2) = round(location(3));
c(1) = round(location(6));
c(2) = round(location(5));
d(1) = round(location(8));
d(2) = round(location(7));
c_p = @(x,y) x(1)*y(2) - x(2)*y(1);
fore2 = zeros(size(foreground));
for xx = 1:x
    for yy = 1:y
        for zz = 1:z
            p = [xx yy];
            ab = a-b;
            bc = b-c;
            cd = c-d;
            da = d-a;
            ap = a-p;
            bp = b-p;
            cp = c-p;
            dp = d-p;
            if (sign(c_p(ab,ap)) > 0 && sign(c_p(bc,bp)) > 0 && sign(c_p(cd,cp)) > 0 && sign(c_p(da,ap)) > 0) ||...
                    (sign(c_p(ab,ap)) < 0 && sign(c_p(bc,bp)) < 0 && sign(c_p(cd,cp)) < 0 && sign(c_p(da,ap)) < 0)
                if foreground(xx,yy,zz) == 1
                    fore2(xx,yy,zz) = 1;
                end
            end
        end
    end
end
im = zeros(x,y,3,z);
for zz = 1:z
    im(:,:,1,zz) = foreground(:,:,zz);
    im(:,:,2,zz) = foreground(:,:,zz) - fore2(:,:,zz);
    im(:,:,3,zz) = foreground(:,:,zz) - fore2(:,:,zz);
end
tifwrite(uint8(im*255),[file_name 'temp_final']);

% remove leakage
load([file_name '\data_c1.mat']);
load([file_name '\data_c2.mat']);
leak_coef = 0.8;
[x, y, z, t] = size(data1);
for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    data1(:,:,:,tt) = data1(:,:,:,tt) - leak_coef*data2(:,:,:,tt);
%     tifwrite(uint8(data1(:,:,:,tt)), [file_name 'temp\' ind]);
end
% data1 = double(uint8(data1));
% save([file_name 'data_c1_processed.mat'], 'data1');

% plot sum of data
% load([file_name '\data_c1_processed.mat']);
% data1 = double(uint8(data1));
data1(data1<0) = 0;
[x, y, z, t] = size(data1);
data = zeros(x,y,z);
for tt = 1:t
    data = data + data1(:,:,:,tt);
end
data = data/t;
save([file_name 'data_c1_processed.mat'], 'data');
tifwrite(uint8(data),[file_name 'temp']);

    


