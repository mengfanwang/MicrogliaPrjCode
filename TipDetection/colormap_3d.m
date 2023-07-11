clc;clear;
% write colormap to 3d tif file

for step = 6:21
step
load(['D:/dropbox/distance_map/' num2str(step) '/1']);

distance_map(distance_map < 0.01) = 0;
distance_map(distance_map>step-0.01) = step-0.01;
distance_map = distance_map/step;
if step == 6
    distance_map_max = distance_map;
    distance_map_arg = ~isnan(distance_map)*6;
else
    distance_map_max = max(distance_map_max,distance_map);
    distance_map_arg(distance_map_max == distance_map) = step;
end
end
distance_map = distance_map_max;
% load prob_31.mat
% distance_map = prob_map;
distance_map = ceil(distance_map*255/max(distance_map(:))-1e-8) + 1;



% distance_map(distance_map>8) = 8;
% distance_map = ceil(distance_map*255/max(distance_map(:))-1e-8) + 1;
color = colormap;
[x,y,z] = size(distance_map);
map = zeros(x,y,3,z); 
for ii = 1:x
    for jj = 1:y
        for kk = 1:z
            if ~isnan(distance_map(ii,jj,kk))
                map(ii,jj,1,kk) = color(distance_map(ii,jj,kk),1);
                map(ii,jj,2,kk) = color(distance_map(ii,jj,kk),2);
                map(ii,jj,3,kk) = color(distance_map(ii,jj,kk),3);
            end
        end
    end
end
tifwrite(uint8(map*255), ['./temp_data/overall']);

                