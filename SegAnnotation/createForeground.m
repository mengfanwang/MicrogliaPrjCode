clc;clear;close all;

folder_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';

for tt = 2:6
    t_ind = num2str(1000+tt);
    t_ind = t_ind(2:end);
    f = tifread([folder_name 'foreground_c2_iter1\' t_ind '.tif']);
    for zz = 1:41
        f_temp = f(:,:,zz);

        ImName = ['D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\gt_seg_initial\' num2str(tt) '_' num2str(zz) '.png'];

        im = zeros(301,301,3);
        im(:,:,1) = f_temp;
        imwrite(im,ImName,'Alpha',f_temp);

    end
end