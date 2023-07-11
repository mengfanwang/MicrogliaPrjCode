% clc;clear;close all;
% file_name = 'D:\dropbox\Modify Series Data-3nd\Priority 1 files\060221 slice1 hippocampus before vessel_Modify Series';

function Registration(file_name)
addpath ./NoRMCorre


load([file_name '\data_c1.mat']);
[x, y, z, t] = size(data1);
options_rigid = NoRMCorreSetParms('d1',x, 'd2',y, 'd3',z, 'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200, 'upd_template', 1);
tic; [M1,shifts1,template1,options_rigid] = normcorre(data1,options_rigid); toc

x_min = floor(min(arrayfun(@(x) x.shifts(1), shifts1)));
x_max = ceil(max(arrayfun(@(x) x.shifts(1), shifts1)));
y_min = floor(min(arrayfun(@(x) x.shifts(2), shifts1)));
y_max = ceil(max(arrayfun(@(x) x.shifts(2), shifts1)));
z_min = floor(min(arrayfun(@(x) x.shifts(3), shifts1)));
z_max = ceil(max(arrayfun(@(x) x.shifts(3), shifts1)));
save([file_name '\drift_c1.mat'], 'shifts1','x_min','x_max','y_min','y_max','z_min','z_max');


%%
% clc;clear;
% file_name = 'D:\dropbox\Modify Series Data\Copy of 032819 slice1 hippo before_Modify Series';
mkdir([file_name '\registration_c1']);
mkdir([file_name '\registration_c2']);
load([file_name '\drift_c1.mat']);
load([file_name '\data_c1.mat']);
load([file_name '\data_c2.mat']);
[x, y, z, t] = size(data1);
coef = 0.8;

x_min = floor(min(arrayfun(@(x) x.shifts(1), shifts1)));
x_max = ceil(max(arrayfun(@(x) x.shifts(1), shifts1)));
y_min = floor(min(arrayfun(@(x) x.shifts(2), shifts1)));
y_max = ceil(max(arrayfun(@(x) x.shifts(2), shifts1)));
z_min = floor(min(arrayfun(@(x) x.shifts(3), shifts1)));
z_max = ceil(max(arrayfun(@(x) x.shifts(3), shifts1)));

data_new1 = zeros(x+x_max-x_min, y+y_max-y_min, z+z_max-z_min,t);
data_new2 = zeros(x+x_max-x_min, y+y_max-y_min, z+z_max-z_min,t);
reg_count = zeros(x+x_max-x_min, y+y_max-y_min, z+z_max-z_min);
for tt = 1:t
    tt
    x_bias = round(shifts1(tt).shifts(1));
    y_bias = round(shifts1(tt).shifts(2));
    z_bias = round(shifts1(tt).shifts(3));
%     data_new1((1:x)+x_bias-x_min,(1:y)+y_bias-y_min,(1:z)+z_bias-z_min,tt) = data1(:,:,:,tt) - coef*data2(:,:,:,tt);
    reg_count((1:x)+x_bias-x_min,(1:y)+y_bias-y_min,(1:z)+z_bias-z_min) = reg_count((1:x)+x_bias-x_min,(1:y)+y_bias-y_min,(1:z)+z_bias-z_min) + 1;
    
    data_new1((1:x)+x_bias-x_min,(1:y)+y_bias-y_min,(1:z)+z_bias-z_min,tt) = data1(:,:,:,tt);
    data_new2((1:x)+x_bias-x_min,(1:y)+y_bias-y_min,(1:z)+z_bias-z_min,tt) = data2(:,:,:,tt);
    
    ind = num2str(1000+tt); 
    ind = ind(2:4);
    tifwrite(uint8(data_new1(:,:,:,tt)),[file_name '\registration_c1\' ind]);
    tifwrite(uint8(data_new2(:,:,:,tt)),[file_name '\registration_c2\' ind]);
    
end
data1 = data_new1;
data2 = data_new2;
save([file_name '\registration_c1'],'data1','reg_count', '-v7.3');
save([file_name '\registration_c2'],'data2', '-v7.3');