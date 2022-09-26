clc;clear;close all;

% read data from czi_file, resize to the same resolution, and save

addpath bfmatlab
file_name = 'D:\dropbox\Modify Series Data\HA-101920-slice1-hippo-vessel-Modify Series';
mkdir(file_name);
mkdir([file_name '\data_c1']);
mkdir([file_name '\data_c2']);
czi_file = bfopen([file_name '.czi']);
czi_file = czi_file{1,1};

% parameter
% czi_file = czi_file(2:2:end,:); % for the second channel
z = 41;
t = 35;
data1 = zeros(301,301,z,t);
data2 = zeros(301,301,z,t);

for tt = 1:t
    for zz = 1:z
        data1(:,:,zz,tt) = cell2mat(czi_file(((tt-1)*z+zz)*2 - 1,1));
        data2(:,:,zz,tt) = cell2mat(czi_file(((tt-1)*z+zz)*2    ,1));
    end
end

% resolution = [0.4613 0.4613 1];
% data1 = interpolation(data1,resolution);
% data2 = interpolation(data2,resolution);


for tt = 1:t
    ind = num2str(1000+tt); 
    ind = [file_name '\data_c1\' ind(2:4)];
    tifwrite(uint8(data1(:,:,:,tt)),ind);
    ind = num2str(1000+tt); 
    ind = [file_name '\data_c2\' ind(2:4)];
    tifwrite(uint8(data2(:,:,:,tt)),ind);
end
save([file_name '\data_c1'],'data1', '-v7.3');
save([file_name '\data_c2'],'data2', '-v7.3');


