clc;clear;close all;
dbstop if error;
% dbstop at 100 if ismember(sub2ind([512 512 61],440,421,1),com_element) 
file_name = 'D:\dropbox\Modify Series Data\Copy of 032619 slice1 hippocampus before_Modify Series';
mkdir([file_name '\registration_c1']);
addpath ./NoRMCorre

% remove leakage
load([file_name '\data_c1.mat']);
[x, y, z, t] = size(data1);

% tic;
% options_rigid = NoRMCorreSetParms('d1',x ,'d2',y, 'd3', z, 'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200, 'upd_template', 1);
% [M1,shifts1,template1,options_rigid] = normcorre(data1,options_rigid);
% toc
% save([file_name '\regist_c1.mat'], 'M1', 'shifts1', 'template1', 'options_rigid');

options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc
save([file_name '\regist_c1_nonrigid.mat'], 'M2', 'shifts2', 'template1', 'options_rigid');
%%
% load([file_name '\regist_c1.mat']);
% M2 = (M1 - min(M1(:))) / (max(M1(:)) - min(M1(:))) * 255;
% for tt = 1:t
%     ind = num2str(1000+tt); 
%     ind = [file_name '\registration_c1\' ind(2:4)];
%     tifwrite(uint8(M2(:,:,:,tt)),ind);
% end