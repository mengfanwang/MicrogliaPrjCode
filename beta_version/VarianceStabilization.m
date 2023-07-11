% clc;clear;close all;
% file_name = 'D:\dropbox\Modify Series Data-3nd\Priority 1 files\060221 slice1 hippocampus before vessel_Modify Series';

function VarianceStabilization(file_name)

channel = 'c2';

mkdir([file_name '\stabilization_' channel]);
load([file_name '\data_' channel]);
if strcmp(channel, 'c1')
    data = data1;
elseif strcmp(channel, 'c2')
    data = data2;
    reg_data = load([file_name '\registration_c2.mat']);
    reg_data = reg_data.data2;
end


kerSize = [3 3 1];
options.histEdges = 0.5:254.5;
options.binEdges = -0.5:255.5;
options.sampleSize = 200;
options.ratio = 0.03;

tic;
[his,variance,parameters] = histogramCount(data,kerSize,options);
toc

[stabilizeFunction, variance] = convexOptimization(his,parameters);

% 0 change to 1 and 255 change to 254
stabilizeFunction(1) = stabilizeFunction(2);
stabilizeFunction(256) = stabilizeFunction(255);
stabilizeFunction = 255*(stabilizeFunction - min(stabilizeFunction))/(max(stabilizeFunction) - min(stabilizeFunction));
save([file_name '\stabilizeFunction_' channel],'stabilizeFunction','variance');

data = reg_data;
data = interp1(0:255,stabilizeFunction,data);
for tt = 1:size(data,4)
    ind = num2str(1000+tt); 
    ind = [file_name '\stabilization_' channel '\' ind(2:4)];
    tifwrite(uint8(data(:,:,:,tt)),ind);
end
save([file_name '\stabilization_' channel '\stabilization_data'], 'data', '-v7.3');