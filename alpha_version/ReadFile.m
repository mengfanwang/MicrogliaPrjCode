clc;clear;close all;

% read data from czi_file, resize to the same resolution, and save
addpath bfmatlab

if isunix
    file_name = '/work/Mengfan/MicrogliaData/100920-slice2-cortex-vessel-after';
else
    file_name = '';
end
mkdir(file_name);
mkdir(fullfile(file_name, 'data_c1'));
mkdir(fullfile(file_name, 'data_c2'));
czi_file = bfopen([file_name '.czi']);
ome = czi_file{1,4};
data_file = czi_file{1,1};

% ablation time detection
x = ome.getPixelsSizeY(0).getValue();
y = ome.getPixelsSizeX(0).getValue();
z = ome.getPixelsSizeZ(0).getValue();
c = ome.getPixelsSizeC(0).getValue();
t = ome.getPixelsSizeT(0).getValue();

ome_data = char(ome.dumpXML());
ome_data = regexp(ome_data,'(?<=DeltaT=")\d+.\d+(?=")','match');
time_second = zeros(length(ome_data)/c/z,1);
for ii = 0:length(ome_data)/c/z-1
    time_second(ii+1) = str2num(ome_data{ii*c*z+1});
end

% save file
data1 = zeros(x,y,z,t);
data2 = zeros(x,y,z,t);

for tt = 1:t
    for zz = 1:z
        data1(:,:,zz,tt) = cell2mat(data_file(((tt-1)*z+zz)*2 - 1,1));
        data2(:,:,zz,tt) = cell2mat(data_file(((tt-1)*z+zz)*2    ,1));
    end
end

ome_data = char(ome.dumpXML());
bit_depth = regexp(ome_data,'(?<=Type="uint)\d+(?=">)','match');
if strcmp(bit_depth{1},'16')
    data1 = data1/257;
    data2 = data2/257;
end
for tt = 1:t
    ind = fullfile(file_name, 'data_c1', num2str(tt));
    tifwrite(uint8(data1(:,:,:,tt)),ind);
    ind = fullfile(file_name, 'data_c2', num2str(tt));
    tifwrite(uint8(data2(:,:,:,tt)),ind);
end
save(fullfile(file_name, 'data_c1'), 'data1', '-v7.3');
save(fullfile(file_name, 'data_c2'), 'data2', '-v7.3');


