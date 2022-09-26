clc;clear;close all;

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';

mkdir([file_name '\stabilization_c2']);
load([file_name '\data_c2']);
data = data2;


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
save([file_name '\stabilizeFunction_c2'],'stabilizeFunction','variance');

data = interp1(0:255,stabilizeFunction,data);
for tt = 1:size(data,4)
    ind = num2str(1000+tt); 
    ind = [file_name '\stabilization_c2\' ind(2:4)];
    tifwrite(uint8(data(:,:,:,tt)),ind);
end
save([file_name '\stabilization_c2\stabilization_data'], 'data');

fiji_descr = ['ImageJ=1.53c' newline ...
            'images=' num2str(size(data,3)*...
                              size(data,4)) newline... 
            'channels=1' newline...
            'slices=' num2str(size(data,3)) newline...
            'frames=' num2str(size(data,4)) newline... 
            'hyperstack=true' newline...
            'mode=grayscale' newline...  
            'loop=false' newline...  
            'min=0.0' newline...      
            'max=255.0']; 
        
t = Tiff([file_name '\stabilization_c2\all.tif'],'w');
tagstruct.ImageLength = size(data,1);
tagstruct.ImageWidth = size(data,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = size(data,3);
tagstruct.Compression = Tiff.Compression.None;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.ImageDescription = fiji_descr;
for tt = 1:size(data,4)
%     for zz = 1:size(data,3)
        t.setTag(tagstruct);
        t.write(uint8(data(:,:,:,tt)));
        t.writeDirectory(); % saves a new page in the tiff file
%     end
end
t.close();

tifwrite(uint8(max(data,[],3)),[file_name '\stabilization_c2\z_prj']);