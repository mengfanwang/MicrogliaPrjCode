clc;clear;
dbstop if error;

file_name = '..\SL-092320-slice1-hippo-vessel-Modify Series';
mkdir([file_name '\foreground']);
load([file_name '\data_c2']);
data = double(data2);
% % 0 change to 1 and 255 change to 254
% stabilizeFunction(1) = stabilizeFunction(2);
% stabilizeFunction(256) = stabilizeFunction(255);
% stabilizeFunction = 255*(stabilizeFunction - min(stabilizeFunction))/(max(stabilizeFunction) - min(stabilizeFunction));
% data = data + 1;
% data = stabilizeFunction(data);
% data = data/sqrt(max(variance));

data_all = data;
foreground = zeros(size(data_all));

for tt = 1:1 %size(data_all,4)
ind = num2str(1000+tt); 
ind = [file_name '\foreground\' ind(2:4)];
data = data_all(:,:,:,tt);


%% parameter setting
smo = 0.8;  % Gaussian parameter
test = imgaussfilt3(randn(200,200,40),smo);
sigma = std(test(:));   % noise estimation
data = imgaussfilt3(data,smo);
dF = data;
thr = 100;    % threshold for region seed
minSize = 50;   % allowed minimum size
Max_rounds = 30; % 30 can be changed. Just an arbitrary number.
minIntensity = 5; % allowed detected minimum intensity

tic;
%% Region seed
labelMap = data>thr;
comp = bwconncomp(labelMap);
comp = comp.PixelIdxList;
sz = cellfun(@numel,comp);

comp = comp(sz>minSize);
arLst = comp;

%% Mask update
labelMap = false(size(dF));
for i = 1:numel(arLst)
    labelMap(arLst{i}) = true;
end

%% Find neighbor pixels

[H,W,T] = size(dF);
candidate_pix = [];
dw = zeros(27,1);
dh = zeros(27,1);
dz = zeros(27,1);
cnt = 1;
for x = -1:1
    for y = -1:1
        for z = -1:1
            dw(cnt) = x;
            dh(cnt) = y;
            dz(cnt) = z;
            cnt = cnt + 1;
        end
    end
end

pix0 = find(labelMap);

for kk = 1:27
    [ih0,iw0,it0] = ind2sub([H,W,T],pix0);
    ih1 = max(min(ih0+dh(kk),H),1);
    iw1 = max(min(iw0+dw(kk),W),1);
    it1 = max(min(it0+dz(kk),T),1);
    pix1 = sub2ind([H,W,T],ih1,iw1,it1);
    candidate_pix = [pix1;candidate_pix];
end
candidate_pix = setdiff(candidate_pix,pix0);
candidate_pix = unique(candidate_pix);

%% Local neighborhood filter
radius = 2;
filter = zeros(radius*2+1,radius*2+1,radius*2+1);
for x = -radius:radius
    for y = -radius:radius
        for z = -radius:radius
            if(x^2 + y^2 + z^2<=radius^2)
                filter(x+radius+1,y+radius+1,z+radius+1) = 1;
            end
        end
    end
end

dF0 = dF;
dF0(~labelMap) = 0;
neighbor_record = imfilter(uint8(labelMap),filter);
summation = imfilter(dF0,filter);
clear dF0;

%% Grow
for nn = 1:Max_rounds
    [pix0,candidate_pix,neighbor_record,summation] = growRegion(pix0,...
        dF,candidate_pix,filter,neighbor_record,summation,sigma,minIntensity);
    fprintf('Grow %d rounds\n',nn);
%     waitbar(nn/30,ff);
end

labelMap(pix0) = true;

comp = bwconncomp(labelMap);
comp = comp.PixelIdxList;
sz = cellfun(@numel,comp);

comp = comp(sz>minSize);
arLst = comp;
% observe just one component
labelMap = false(size(dF));
for i = 1:numel(arLst)
    labelMap(arLst{i}) = true;
end
toc;

%% Output

% arLst = bwconncomp(labelMap);
% arLst = arLst.PixelIdxList;
% ov1 = plt.regionMapWithData(arLst,0.5*dF/max(dF(:)),0.5,[]);
% ov = uint8(zeros(H,W,T));
% ov = uint8(max(labelMap*255,data*sqrt(max(variance))));
ov = uint8(labelMap*255);
% ov(:,:,2,:) = uint8(labelMap*255);
% ov(:,:,3,:) = uint8(labelMap*255);
% ov(:,1:W,1,:) = uint8(dF/max(dF(:))*255);
% ov(:,1:W,2,:) = uint8(dF/max(dF(:))*255);
% ov(:,1:W,3,:) = uint8(dF/max(dF(:))*255);
% dF0 = dF/max(dF(:));
% dF0 = max(dF0,double(labelMap));
foreground(:,:,:,tt) = labelMap;
tifwrite(ov,ind);
end
foreground = logical(foreground);
save([file_name '\foreground'],'foreground');

