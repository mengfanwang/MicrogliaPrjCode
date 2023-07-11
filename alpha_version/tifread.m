function ImStack = tifread(impath)
% modified by Mengfan
% read a tif file
% impath: tif path e.g. 'C:\\a.tif'
% can read the whole image only
% support both grayscale and RGB 3Dimage
info = imfinfo(impath);
num_images = numel(info);
org_h = info(1).Height;
org_w = info(1).Width;
channel = info(1).PhotometricInterpretation;

if strcmp(channel, 'BlackIsZero')
    ImStack = zeros(org_h,org_w,num_images);
    for i = 1:num_images
        Im = imread(impath,i);
        ImStack(:,:,i) = double(Im);    
    end
elseif strcmp(channel, 'RGB')
    ImStack = zeros(org_h,org_w,3,num_images);
    for i = 1:num_images
        Im = imread(impath,i);
        ImStack(:,:,:,i) = double(Im);    
    end
end
    