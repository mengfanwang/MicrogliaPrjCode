clc;clear;close all;
dbstop if error

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';
cur_file = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series.czi';
addpath ./bfmatlab
load([file_name 'data_c2.mat']);
vid = data2;
vesselFlag = 0;

[h,w, z, t] = size(vid);
%% data resolution
abSite.resolution = [];
data = bfopen(cur_file);
ome = data{1,4};
% ome = loci.formats.MetadataTools.createOMEXMLMetadata();
% reader.setMetadataStore(ome);
% reader.setId(cur_file);
rx = double(ome.getPixelsPhysicalSizeX(0).value);
ry = double(ome.getPixelsPhysicalSizeY(0).value);
rz = double(ome.getPixelsPhysicalSizeZ(0).value);
totaltime = double(ome.getPlaneDeltaT(0,z*t-1).value) - double(ome.getPlaneDeltaT(0,0).value);
abSite.frameAvgInterval = (totaltime / (t-1))/60;
abSite.frameIntervals = zeros(t-1,1);
if vesselFlag
    for i=1:t-1
        abSite.frameIntervals(i) = ...
            double(ome.getPlaneDeltaT(0,2*i*z).value) - double(ome.getPlaneDeltaT(0,2*(i-1)*z).value);
    end
else
    for i=1:t-1
        abSite.frameIntervals(i) = ...
            double(ome.getPlaneDeltaT(0,i*z).value) - double(ome.getPlaneDeltaT(0,(i-1)*z).value);
    end
end
abSite.frameIntervals = abSite.frameIntervals./60;
% reader.close();
abSite.resolution = [rx, ry, rz];

%% ablation time
% maxProj = max(vid,[], 3);
% cv = zeros(1,t);
% for i=1:t
%     crop = maxProj(max(1,round(locy)-2):min(h,round(locy)+2), ...
%         max(1,round(locx)-2):min(w,round(locx)+2),i);
%     cv(i) = nanmean(crop(:));
% end
% if vesselFlag
%     [~, pos] = max(cv(2:end)-cv(1:end-1));
% else
%     [~, pos] = min(cv(2:end)-cv(1:end-1));
% end
[~, pos] = max(abSite.frameIntervals);
abSite.timePoint = pos;

%% ablation site
% x and y location
zz = czifinfo( cur_file ); 
loc0 = strfind(zz.metadataXML, '<CenterX>');
loc1 = strfind(zz.metadataXML, '</CenterX>');
if ~isempty(loc0) && ~isempty(loc1)
    locx = str2double(zz.metadataXML(loc0(1)+9:loc1(1)-1));

    loc0 = strfind(zz.metadataXML, '<CenterY>');
    loc1 = strfind(zz.metadataXML, '</CenterY>');
    locy = str2double(zz.metadataXML(loc0(1)+9:loc1(1)-1));
else
    %
    disp('No center captured, estimate needed');
    loc0 = strfind(zz.metadataXML, '<Points>');
    loc1 = strfind(zz.metadataXML, '</Points>');
    strTmp = zz.metadataXML(loc0(1)+8:loc1(1)-1);
    loc2 = strfind(strTmp, ' ');
    loc3 = strfind(strTmp, ',');
    
    locx = str2double(strTmp(loc2(1)+1:loc3(2)-1));
    locy = str2double(strTmp(loc3(2)+1:loc2(2)-1));

end

% z-location
if vesselFlag
    testVid = vid(:,:,:,pos+1);
else
    testVid = vid(:,:,:,pos);
end
cvZ = zeros(size(testVid,3),1);
for i=1:size(testVid,3)
    crop = testVid(max(1,round(locy)-2):min(h,round(locy)+2), ...
        max(1,round(locx)-2):min(w,round(locx)+2),i);
    cvZ(i) = mean(crop(:));
end
[~, locz] = max(cvZ);
abSite.loc = [locx, locy, locz];