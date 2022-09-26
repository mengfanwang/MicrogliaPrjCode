clc;clear;close all;

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
mkdir([file_name '\soma']);
load([file_name '\data_c2.mat']);
% tic;
% for tt = 1:35
%     tt
%     ind = num2str(1000+tt);
%     ind = ind(2:4);
% %     bw = imbinarize(data2(:,:,:,tt));
%     bw = data2(:,:,:,tt) > 100;
%     tifwrite(uint8(bw*255),[file_name '\soma\' ind]);
% end
% toc


for tt = 1:35
    tt
ind = num2str(1000+tt);
ind = ind(2:4);
% data = load([file_name '\foreground_c2_reconnect\' ind '.mat']);
% data = data.foreground;
data = data2(:,:,:,tt) > 100;
[x,y,z] = size(data); 

% use ellipsoid to find soma
data_radius = zeros(x,y,z);
threshold_intensity = 0.6;  % the average intensity in ellipsoid
continue_flag = 1;
radius = 1; % radius is the z-radius
while continue_flag
    % build kernel
    kernel = build_kernel(radius);
    data_radius_temp = convn(data,kernel,'same');
    data_radius_temp = data_radius_temp/sum(kernel(:));
    data_radius_temp = data_radius_temp > threshold_intensity;
    if sum(data_radius_temp(:)) == 0
        continue_flag = 0;
    else
        data_radius(data_radius_temp) = radius;
        radius = radius + 1;
    end
end

% remove small components
% load('data_radius.mat');
threshold_radius = 3;
components = bwconncomp(data_radius,26);
components_len = length(components.PixelIdxList);
final_result = zeros(x,y,z);


for com_ind = 1:components_len
    locations = components.PixelIdxList{1,com_ind};
    [max_radius,max_ind] = max(data_radius(locations));
    [x_center, y_center, z_center] = ind2sub([x,y,z],locations(max_ind));  
    if max_radius > threshold_radius
        if length(locations) > 100
            for ii = max(1,x_center-2*max_radius):min(x,x_center+2*max_radius)
                for jj = max(1,y_center-2*max_radius):min(y,y_center+2*max_radius)
                    for kk = max(1,z_center-max_radius):min(z,z_center+max_radius)
                        if (ii-x_center)^2/(2*max_radius)^2+(jj-y_center)^2/(2*max_radius)^2+(kk-z_center)^2/max_radius^2 <= 1 && data(ii,jj,kk) == 1
                            final_result(ii,jj,kk) = 255;
                        end
                    end
                end
            end
        end
    end
end
final_result = uint8(final_result);


tifwrite(uint8(final_result*255),[file_name '\soma\' ind]);
soma = final_result;
soma = logical(soma);
save([file_name '\soma\' ind '.mat'],'soma');
end


function kernel = build_kernel(radius)
    kernel = zeros(4*radius+1,4*radius+1,2*radius+1);
    for ii = 1:4*radius+1
        for jj = 1:4*radius+1
            for kk = 1:2*radius+1
                if (ii-2*radius-1)^2/(2*radius)^2+(jj-2*radius-1)^2/(2*radius)^2+(kk-radius-1)^2/radius^2 <= 1
                    kernel(ii,jj,kk) = 1;
                end
            end
        end
    end
end




