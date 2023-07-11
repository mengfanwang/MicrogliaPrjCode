% clc;clear;close all;
% dbstop if error
% 
% file_name = 'D:\dropbox\Modify Series Data\Copy of 042519 slice1 hippo before_Modify Series';

function Output(file_name, path)
addpath ./bfmatlab

file_split = split(file_name,'\');
if ~isfolder([path 'output_new'])
    mkdir([path 'output_new']);
end
xls_name = [path 'output_new\' file_split{end} '.xlsx'];


% xls_name = 'Copy of 042519 slice1 hippo before_Modify Series.xlsx';
xlswrite(xls_name, {'Path'}, 'Sheet1','A1:A1');
xlswrite(xls_name, {file_split{end}}, 'Sheet1','B1:B1');

% basic information
czi_file = bfopen([file_name '.czi']);
ome = czi_file{1,4};
% rx = double(ome.getPixelsPhysicalSizeX(0).value);
% ry = double(ome.getPixelsPhysicalSizeY(0).value);
% rz = double(ome.getPixelsPhysicalSizeZ(0).value);
rx = 0.46125879; ry = rx; rz = 1;
unit_volume = rx*ry*rz;
dv = [rx ry rz];
x = czi_file{1,4}.getPixelsSizeX(0).getValue();
y = czi_file{1,4}.getPixelsSizeY(0).getValue();
z = czi_file{1,4}.getPixelsSizeZ(0).getValue();
c = czi_file{1,4}.getPixelsSizeC(0).getValue();
% t = czi_file{1,4}.getPixelsSizeT(0).getValue();
total_volume = x*y*z*unit_volume;

% ablation time detection
ome_data = char(czi_file{1,4}.dumpXML());
ome_data = regexp(ome_data,'(?<=DeltaT=")\d+.\d+(?=")','match');
time_second = zeros(length(ome_data)/c/z,1);
for ii = 0:length(ome_data)/c/z-1
    time_second(ii+1) = str2num(ome_data{ii*c*z+1});
end
[~, ablation_time] = max(diff(time_second));
ablation_time = ablation_time + 1;

xlswrite(xls_name, {'Ablation time'},'Sheet1','A2:A2');
xlswrite(xls_name, time_second(ablation_time),'Sheet1','B2:B2');
xlswrite(xls_name, {'Ablation frame'},'Sheet1','A3:A3');
xlswrite(xls_name, ablation_time,'Sheet1','B3:B3');

% ablation region detection
% location = [140.911357 158.143121 139.799631 163.701754 160.088643 168.426593 161.756233 160.088643];
load([file_name '\drift_c1.mat']);
load([file_name '\foreground_c1.mat']);
foreground_c1 = foreground;
[x,y,z] = size(foreground);

% % Prority 1/2/3/5 use polygons
% location = getAblationPoint(file_name);
% location(:,1) = location(:,1) + shifts1(1).shifts(1) - x_min;
% location(:,2) = location(:,2) + shifts1(1).shifts(2) - y_min;
% ablation_region = getAblationRegion(foreground_c1, location);
% [~,~,z_ind] = ind2sub(size(ablation_region),find(ablation_region));
% ablation_center = mean(location);
% ablation_center = [ablation_center mean(z_ind)];

% Prority 4 use circles
location = getAblationCircle(file_name);
location(1) = location(1) + shifts1(1).shifts(1) - x_min;
location(2) = location(2) + shifts1(1).shifts(2) - y_min;
ind = num2str(1000+ablation_time); 
ind = ind(2:4);
load([file_name '.\soma\' ind '.mat'],'soma');
ablation_region = getAblationRegion_Circle(soma, location);
[~,~,z_ind] = ind2sub(size(ablation_region),find(ablation_region));
if isempty(z_ind)
    z_ind = z/2 + shifts1(1).shifts(3) - z_min;
end
ablation_center = location(1:2);
ablation_center = [ablation_center mean(z_ind)];

im_vessel = zeros(x,y,3,z);
for zz = 1:z
    im_vessel(:,:,1,zz) = foreground_c1(:,:,zz);
    im_vessel(:,:,2,zz) = foreground_c1(:,:,zz) - ablation_region(:,:,zz);
    im_vessel(:,:,3,zz) = foreground_c1(:,:,zz) - ablation_region(:,:,zz);
%     im_vessel(:,:,1,zz) = ablation_region(:,:,zz);
end
tifwrite(uint8(im_vessel*255),[file_name '\ablation_site']);

if exist([file_name '\distance_c1.mat'],'file')
    load([file_name '\distance_c1']);
else
    tic;
%     if sum(foreground_c1(:)) > 0
        distance_c1 = getFrameDist(foreground_c1, dv);
%     else
%         distance_c1 = zeros(size(foreground_c1));
%     end
    toc
    save([file_name '\distance_c1'],'distance_c1');   
end
 

% find 2 min and 10 min
time_minute = (time_second - time_second(ablation_time))/60;
[~, frame_2min] = min(abs(time_minute-2));
[~, frame_6min] = min(abs(time_minute-6));
[~, frame_10min] = min(abs(time_minute-10));
t = min(frame_10min + 2,czi_file{1,4}.getPixelsSizeT(0).getValue());
xlswrite(xls_name,{'Frame after 2 min'},'Sheet1','A4:A4');
xlswrite(xls_name,frame_2min,'Sheet1','B4:B4');
xlswrite(xls_name,{'Frame after 6 min'},'Sheet1','A5:A5');
xlswrite(xls_name,frame_6min,'Sheet1','B5:B5');
xlswrite(xls_name,{'Frame after 10 min'},'Sheet1','A6:A6');
xlswrite(xls_name,frame_10min,'Sheet1','B6:B6');
foreground_c1_num = sum(foreground_c1(:));

xlswrite(xls_name,{'Vessel volume'},'Sheet1','A9:A9');
xlswrite(xls_name,{'Glia volume(begining)'},'Sheet1','A10:A10');
xlswrite(xls_name,{'Glia volume(ablation)'},'Sheet1','A11:A11');
xlswrite(xls_name,{'Glia volume(2min)'},'Sheet1','A12:A12');
xlswrite(xls_name,{'Glia volume(6min)'},'Sheet1','A13:A13');
xlswrite(xls_name,{'Glia volume(10min)'},'Sheet1','A14:A14');
xlswrite(xls_name,{'Glia volume(end)'},'Sheet1','A15:A15');
xlswrite(xls_name,foreground_c1_num*unit_volume,'Sheet1','B9:B9');

frame_list = [1 ablation_time frame_2min frame_6min frame_10min t];
for ii = 1:6
    tt = frame_list(ii);
    ind = num2str(1000+tt);
    ind = ind(2:4);
%     load([file_name '\foreground_c2_reconnect\' ind '.mat']); 
    load([file_name '\foreground_c2\' ind '.mat']);
    foreground_c2 = foreground;
    if ii == 2
        foreground_ablation = foreground;
    end
    foreground_c2_num = sum(foreground_c2(:));
    
    jj = num2str(9+ii);
    xlswrite(xls_name,foreground_c2_num*unit_volume,'Sheet1',['B' jj ':B' jj]);
end
xlswrite(xls_name,{'Total volume of view'},'Sheet1','A16:A16');
xlswrite(xls_name,total_volume,'Sheet1','B16:B16');


%%% soma
load([file_name '\stabilization_c2\stabilization_data.mat']);
data_all = data;
[x, y, z, t] = size(data_all);
soma_all = zeros(x,y,z);
for tt = 1:t
    ind = num2str(1000+tt); 
    ind = ind(2:4);
    load([file_name '.\soma\' ind '.mat'],'soma');
    soma_all = soma_all + soma;
end
soma_all = soma_all > 0;
soma_all = imerode(soma_all,ones(3,3,3));
tifwrite(uint8(soma_all*255),[file_name '\soma']);
soma_components = bwconncomp(soma_all);
soma_num = 0;
xlswrite(xls_name,{'Num of cell bodies at ablation'},'Sheet1','A19:A19');
xlswrite(xls_name,{'Cell body ID'},'Sheet1','A20:A20');
xlswrite(xls_name,{'distance2vessel'},'Sheet1','B20:B20');
xlswrite(xls_name,{'distance2ablation_site'},'Sheet1','C20:C20');
for ss = 1:soma_components.NumObjects
%     if length(soma_components.PixelIdxList{ss}) > 50
    soma_num = soma_num + 1;
    [x_ind, y_ind, z_ind] = ind2sub(size(soma), soma_components.PixelIdxList{ss});
    x_soma = mean(x_ind); y_soma = mean(y_ind); z_soma = mean(z_ind);
    soma2vessel = distance_c1(round(x_soma),round(y_soma),round(z_soma));
    soma2ablation = getP2Vdist([x_soma y_soma z_soma],ablation_center,dv);
    jj = num2str(20+soma_num);
    xlswrite(xls_name,{num2str(ss)},'Sheet1',['A' jj ':A' jj]);
    xlswrite(xls_name,soma2vessel,'Sheet1',['B' jj ':B' jj]);
    xlswrite(xls_name,soma2ablation,'Sheet1',['C' jj ':C' jj]);
%     end
end
xlswrite(xls_name,soma_num,'Sheet1','B19:B19');


%% tracking
load([file_name '\tip_info.mat']);
addpath ./CINDA
v_max = 10;
num_tip = cellfun(@length,xCoord);
% cum_tip = cumsum(num_tip);
detection_arcs = zeros(sum(num_tip),4);
node_table = zeros(sum(num_tip),4);
node_matrix = zeros(size(data_all));
transition_arcs = zeros(sum(num_tip(1:t-1).*num_tip(2:t)),3);
transition_ind = 0;
for tt = 1:t
    for jj = 1:num_tip(tt)
        ind = sum(num_tip(1:tt-1))+jj;
        detection_arcs(ind,:) = [ind v_max/2 v_max/2 -v_max];
        node_table(ind,:) = [tt xCoord{tt}(jj) yCoord{tt}(jj) zCoord{tt}(jj)];
        node_matrix(xCoord{tt}(jj),yCoord{tt}(jj),zCoord{tt}(jj),tt) = ind;
        if tt ~= t
            for kk = 1:num_tip(tt+1)
                transition_ind = transition_ind + 1;
                transition_arcs(transition_ind,:) = [ind sum(num_tip(1:tt))+kk ...
                    sqrt((xCoord{tt}(jj)-xCoord{tt+1}(kk))^2 + (yCoord{tt}(jj)-yCoord{tt+1}(kk))^2 + (zCoord{tt}(jj)-zCoord{tt+1}(kk))^2)];
            end
        end
    end
end

[trajectories, costs] = mcc4mot(detection_arcs,transition_arcs);
trajectories(cellfun(@length,trajectories)<4) = [];

%% Sheet 2 new variable
% xlswrite([file_name '\' xls_name],{'track_ID'},'Sheet2','A1:A1');
% xlswrite([file_name '\' xls_name],{'first_frame'},'Sheet2','B1:B1');
% xlswrite([file_name '\' xls_name],{'last_frame'},'Sheet2','C1:C1');
% xlswrite([file_name '\' xls_name],{'total_distance'},'Sheet2','D1:D1');
% xlswrite([file_name '\' xls_name],{'end_distance'},'Sheet2','E1:E1');
% xlswrite([file_name '\' xls_name],{'tortuosity'},'Sheet2','F1:F1');
% xlswrite([file_name '\' xls_name],{'valid_at_ablation'},'Sheet2','G1:G1');
% xlswrite([file_name '\' xls_name],{'distance2vessel_at_ablation'},'Sheet2','H1:H1');
% xlswrite([file_name '\' xls_name],{'cell_body_identity'},'Sheet2','I1:I1');
track_ID = zeros(length(trajectories),1);
first_frame = zeros(length(trajectories),1);
last_frame = zeros(length(trajectories),1);
total_distance = zeros(length(trajectories),1);
end_distance = zeros(length(trajectories),1);
tortuosity = zeros(length(trajectories),1);
valid_at_ablation = zeros(length(trajectories),1);
distance2vessel_at_ablation = zeros(length(trajectories),1);
cell_body_identity = zeros(length(trajectories),1);
cell_body_temp = zeros(length(trajectories),soma_components.NumObjects);
% label_ablation = bwlabeln(foreground_ablation, 26);
% load foreground2soma distance for every frame and every soma
tic;
distance_soma = zeros(x,y,z,t,soma_components.NumObjects);
for tt = 1:t
    tt
    ind = num2str(1000+tt);
    ind = ind(2:4);
%     load([file_name '\foreground_c2_reconnect\' ind '.mat']);
    load([file_name '\foreground_c2\' ind '.mat']);
    foreground_c2 = foreground;
    for ii = 1:soma_components.NumObjects
        soma_mask = zeros(x,y,z);
        soma_mask(soma_components.PixelIdxList{ii}) = 1;
        distance_soma(:,:,:,tt,ii) = bwdistgeodesic(logical(foreground_c2),logical(soma_mask),'quasi-euclidean');
    end
end
toc

for ii = 1:length(trajectories)
    ii
    track_ID(ii) = ii;
    start_ind = trajectories{ii}(1);
    end_ind = trajectories{ii}(end);
    first_frame(ii) = node_table(start_ind,1);
    last_frame(ii) = node_table(end_ind,1);
    
%     total_distance = 0;
    end_distance(ii) = getP2Pdist(node_table(end_ind,2:4), node_table(start_ind,2:4), dv);
    for jj = 2:length(trajectories{ii})
        total_distance(ii) = total_distance(ii) + getP2Pdist(node_table(trajectories{ii}(jj-1),2:4), node_table(trajectories{ii}(jj),2:4), dv);
    end
    tortuosity(ii) = total_distance(ii)/end_distance(ii);
    
    t_temp = zeros(length(trajectories{ii}),1);
    for jj = 1:length(trajectories{ii})
        t_temp(jj) = node_table(trajectories{ii}(jj),1);
    end
    if ismember(ablation_time, t_temp)
        valid_at_ablation(ii) = 1;
        [~, frame_loc] = ismember(ablation_time, t_temp);
        frame_ind = trajectories{ii}(frame_loc);
        distance2vessel_at_ablation(ii) = distance_c1(node_table(frame_ind,2),node_table(frame_ind,3),node_table(frame_ind,4));
    else
        valid_at_ablation(ii) = 0;
        distance2vessel_at_ablation(ii) = -1;
    end
    cell_body_identity(ii) = 0;
    for jj = 1:length(trajectories{ii})
        frame_ind = trajectories{ii}(jj);
        cell_body_temp_dist = zeros(1,soma_components.NumObjects);
        for kk = 1:soma_components.NumObjects
            cell_body_temp_dist(kk) = distance_soma(node_table(frame_ind,2),node_table(frame_ind,3),node_table(frame_ind,4),node_table(frame_ind,1),kk);            
        end
        [~,kk] = min(cell_body_temp_dist);
        cell_body_temp(ii,kk) = cell_body_temp(ii,kk) + 1;
    end
    [~, soma_max] = max(cell_body_temp(ii,:));
    if cell_body_temp(ii,soma_max) > length(trajectories{ii})/2
        cell_body_identity(ii) = soma_max;
    end
%     jj = num2str(ii+1);
%     xlswrite([file_name '\' xls_name],ii,'Sheet2',['A' jj ':A' jj]);
%     xlswrite([file_name '\' xls_name],first_frame,'Sheet2',['B' jj ':B' jj]);
%     xlswrite([file_name '\' xls_name],last_frame,'Sheet2',['C' jj ':C' jj]);
%     xlswrite([file_name '\' xls_name],total_distance,'Sheet2',['D' jj ':D' jj]);
%     xlswrite([file_name '\' xls_name],end_distance,'Sheet2',['E' jj ':E' jj]);
%     xlswrite([file_name '\' xls_name],total_distance/end_distance,'Sheet2',['F' jj ':F' jj]);
%     xlswrite([file_name '\' xls_name],valid_at_ablation,'Sheet2',['G' jj ':G' jj]);
%     xlswrite([file_name '\' xls_name],distance2vessel_at_ablation,'Sheet2',['H' jj ':H' jj]);
%     xlswrite([file_name '\' xls_name],cell_body_identity,'Sheet2',['I' jj ':I' jj]);
end
table2 = table(track_ID,first_frame,last_frame,total_distance,end_distance,...
    tortuosity,valid_at_ablation,distance2vessel_at_ablation,cell_body_identity);
writetable(table2,xls_name,'Sheet',2);

%% Sheet 3 old variable
% xlswrite([file_name '\' xls_name],{'track_ID'},'Sheet3','A1:A1');
% xlswrite([file_name '\' xls_name],{'distance_baseline'},'Sheet3','B1:B1');
% xlswrite([file_name '\' xls_name],{'speed_baseline'},'Sheet3','C1:C1');
% xlswrite([file_name '\' xls_name],{'distance_2min_after_ablation'},'Sheet3','D1:D1');
% xlswrite([file_name '\' xls_name],{'speed_2min_after_ablation'},'Sheet3','E1:E1');
% xlswrite([file_name '\' xls_name],{'distance_6min_after_ablation'},'Sheet3','F1:F1');
% xlswrite([file_name '\' xls_name],{'speed_6min_after_ablation'},'Sheet3','G1:G1');
% xlswrite([file_name '\' xls_name],{'distance_10min_after_ablation'},'Sheet3','H1:H1');
% xlswrite([file_name '\' xls_name],{'speed_10min_after_ablation'},'Sheet3','I1:I1');
% xlswrite([file_name '\' xls_name],{'distance2ablation_start'},'Sheet3','J1:J1');
% xlswrite([file_name '\' xls_name],{'distance2ablation_at_ablation'},'Sheet3','K1:K1');
% xlswrite([file_name '\' xls_name],{'distance2ablation_end'},'Sheet3','L1:L1');
track_ID = zeros(length(trajectories),1);
distance_baseline = zeros(length(trajectories),1);
speed_baseline = zeros(length(trajectories),1);
distance_2min_after_ablation = zeros(length(trajectories),1);
speed_2min_after_ablation = zeros(length(trajectories),1);
distance_6min_after_ablation = zeros(length(trajectories),1);
speed_6min_after_ablation = zeros(length(trajectories),1);
distance_10min_after_ablation = zeros(length(trajectories),1);
speed_10min_after_ablation = zeros(length(trajectories),1);
distance2ablation_start = zeros(length(trajectories),1);
distance2ablation_at_ablation = zeros(length(trajectories),1);
distance2ablation_end = zeros(length(trajectories),1);
for ii = 1:length(trajectories)
    ii
    track_ID(ii) = ii;
    t_temp = zeros(length(trajectories{ii}),1);
    for jj = 1:length(trajectories{ii})
        t_temp(jj) = node_table(trajectories{ii}(jj),1);
    end
    start_ind = trajectories{ii}(1);
    end_ind = trajectories{ii}(end);
    if (t_temp(1) < ablation_time) && ismember(ablation_time, t_temp)
        distance_baseline(ii) = getTrackDist(trajectories{ii}, node_table, t_temp(1), ablation_time, dv);
        speed_baseline(ii) = distance_baseline(ii)/(time_second(ablation_time) - time_second(t_temp(1))) * 60;
    else
        distance_baseline(ii) = -1;
        speed_baseline(ii) = -1;
    end
    
    if ismember(ablation_time, t_temp) && ismember(frame_2min, t_temp)
        distance_2min_after_ablation(ii) = getTrackDist(trajectories{ii}, node_table, ablation_time, frame_2min, dv);
    else
        distance_2min_after_ablation(ii) = -1;
    end
    if ismember(frame_2min, t_temp)
        speed_2min_after_ablation(ii) = getSpeed(trajectories{ii}, node_table, frame_2min, dv, time_second) * 60;
    else
        speed_2min_after_ablation(ii) = -1;
    end
    
    if ismember(ablation_time, t_temp) && ismember(frame_6min, t_temp)
        distance_6min_after_ablation(ii) = getTrackDist(trajectories{ii}, node_table, ablation_time, frame_6min, dv);
    else
        distance_6min_after_ablation(ii) = -1;
    end
    if ismember(frame_6min, t_temp)
        speed_6min_after_ablation(ii) = getSpeed(trajectories{ii}, node_table, frame_6min, dv, time_second) * 60;
    else
        speed_6min_after_ablation(ii) = -1;
    end
    
    if ismember(ablation_time, t_temp) && ismember(frame_10min, t_temp)
        distance_10min_after_ablation(ii) = getTrackDist(trajectories{ii}, node_table, ablation_time, frame_10min, dv);
    else
        distance_10min_after_ablation(ii) = -1;
    end
    if ismember(frame_10min, t_temp)
        speed_10min_after_ablation(ii) = getSpeed(trajectories{ii}, node_table, frame_10min, dv, time_second) * 60;
    else
        speed_10min_after_ablation(ii) = -1;
    end
    
    distance2ablation_start(ii) = getP2Pdist(node_table(start_ind,2:4), ablation_center, dv);
    distance2ablation_end(ii) = getP2Pdist(node_table(end_ind,2:4), ablation_center, dv);    
    if ismember(ablation_time, t_temp)
        [~, frame_loc] = ismember(ablation_time, t_temp);
        frame_ind = trajectories{ii}(frame_loc);
        distance2ablation_at_ablation(ii) = getP2Pdist(node_table(frame_ind,2:4), ablation_center, dv);  
    else
        distance2ablation_at_ablation(ii) = -1;
    end
    
    
%     jj = num2str(ii+1);
%     xlswrite([file_name '\' xls_name],ii,'Sheet3',['A' jj ':A' jj]);
%     xlswrite([file_name '\' xls_name],distance_baseline,'Sheet3',['B' jj ':B' jj]);
%     xlswrite([file_name '\' xls_name],speed_baseline,'Sheet3',['C' jj ':C' jj]);
%     xlswrite([file_name '\' xls_name],distance_2min,'Sheet3',['D' jj ':D' jj]);
%     xlswrite([file_name '\' xls_name],speed_2min,'Sheet3',['E' jj ':E' jj]);
%     xlswrite([file_name '\' xls_name],distance_6min,'Sheet3',['F' jj ':F' jj]);
%     xlswrite([file_name '\' xls_name],speed_6min,'Sheet3',['G' jj ':G' jj]);
%     xlswrite([file_name '\' xls_name],distance_10min,'Sheet3',['H' jj ':H' jj]);
%     xlswrite([file_name '\' xls_name],speed_10min,'Sheet3',['I' jj ':I' jj]);
%     xlswrite([file_name '\' xls_name],distance2ablation_start,'Sheet3',['J' jj ':J' jj]);
%     xlswrite([file_name '\' xls_name],distance2ablation_ablation,'Sheet3',['K' jj ':K' jj]);
%     xlswrite([file_name '\' xls_name],distance2ablation_end,'Sheet3',['L' jj ':L' jj]);
end
table3 = table(track_ID, distance_baseline, speed_baseline, distance_2min_after_ablation, speed_2min_after_ablation,...
    distance_6min_after_ablation, speed_6min_after_ablation, distance_10min_after_ablation, speed_10min_after_ablation,...
    distance2ablation_start, distance2ablation_at_ablation, distance2ablation_end);
writetable(table3,xls_name,'Sheet',3);

end

function distance = getFrameDist(foreground, dv)
    [x_ind, y_ind, z_ind] = ind2sub(size(foreground), find(foreground));
    [x, y, z] = size(foreground);
    distance = zeros(size(foreground));
    for xx = 1:x
        xx
        for yy = 1:y
            for zz = 1:z
                distance(xx,yy,zz) = min(getP2Vdist([xx yy zz], [x_ind y_ind z_ind], dv));
            end
        end
    end
end

function distance = getP2Vdist(p1,p2,dv)
    distance = zeros(size(p2,1),1);
    for ii = 1:3
        distance = distance + ((p1(ii)-p2(:,ii)).*dv(ii)).^2;
    end
    distance = sqrt(distance);
end

function distance = getP2Pdist(p1,p2,dv)
    distance = 0;
    for ii = 1:3
        distance = distance + ((p1(ii)-p2(ii))*dv(ii))^2;
    end
    distance = sqrt(distance);
end

function speed = getSpeed(trajectory, node_table, frame, dv, time_second)
    
    t_temp = zeros(length(trajectory),1);
    for jj = 1:length(trajectory)
        t_temp(jj) = node_table(trajectory(jj),1);
    end
    start_ind = trajectory(1);
    end_ind = trajectory(end);
    if node_table(end_ind,1) < frame || node_table(start_ind,1) > frame
        speed = -1;
    elseif node_table(end_ind,1) == frame
        start_ind = trajectory(end-1);
        distance = getP2Pdist(node_table(end_ind,2:4), node_table(start_ind,2:4), dv);
        speed = distance/(time_second(node_table(end_ind,1))-time_second(node_table(start_ind,1)));
    elseif node_table(start_ind,1) == frame
        end_ind = trajectory(2);
        distance = getP2Pdist(node_table(end_ind,2:4), node_table(start_ind,2:4), dv);
        speed = distance/(time_second(node_table(end_ind,1))-time_second(node_table(start_ind,1)));
    elseif ismember(frame, t_temp)
        [~, frame_loc] = ismember(frame, t_temp);
        start_ind = trajectory(frame_loc-1);
        end_ind = trajectory(frame_loc);
        distance = getP2Pdist(node_table(end_ind,2:4), node_table(start_ind,2:4), dv);
        start_ind = trajectory(frame_loc);
        end_ind = trajectory(frame_loc+1);
        distance = distance + getP2Pdist(node_table(end_ind,2:4), node_table(start_ind,2:4), dv);
        speed = distance/(time_second(node_table(trajectory(frame_loc+1),1))-time_second(node_table(trajectory(frame_loc-1),1)));
    else
        error('frame not find');
    end
end

function distance = getTrackDist(trajectory, node_table, frame1, frame2, dv)
    distance = 0;
    for jj = 1:length(trajectory)
        t_temp = node_table(trajectory(jj),1);
        if t_temp >= frame1 && t_temp < frame2
            start_ind = trajectory(jj);
            end_ind = trajectory(jj+1);
            distance = distance + getP2Pdist(node_table(end_ind,2:4), node_table(start_ind,2:4), dv);
        end
    end
end

function location = getAblationPoint(file_name)
    zz = czifinfo([file_name '.czi']); 
    loc0 = strfind(zz.metadataXML, '<Points>');
    loc1 = strfind(zz.metadataXML, '</Points>');
    strTmp = zz.metadataXML(loc0(1)+8:loc1(1)-1);
    loc2 = strfind(strTmp, ' ');
    loc3 = strfind(strTmp, ',');
    num_points = length(loc3);
    location = zeros(num_points,2);
    location(1,2) = str2double(strTmp(1:loc3(1)-1));
    location(1,1) = str2double(strTmp(loc3(1)+1:loc2(1)-1)); 
    for ii = 2:num_points - 1
        location(ii,2) = str2double(strTmp(loc2(ii-1)+1:loc3(ii)-1));
        location(ii,1) = str2double(strTmp(loc3(ii)+1:loc2(ii)-1));
    end
    location(num_points,2) = str2double(strTmp(loc2(num_points-1)+1:loc3(num_points)-1));
    location(num_points,1) = str2double(strTmp(loc3(num_points)+1:end));
end

function location = getAblationCircle(file_name)
    % locataion: [x y radius]
    zz = czifinfo([file_name '.czi']); 
    location = zeros(1,3);
    
    loc0 = strfind(zz.metadataXML, '<CenterX>');
    loc1 = strfind(zz.metadataXML, '</CenterX>');
    strTmp = zz.metadataXML(loc0(1)+9:loc1(1)-1);
    location(2) = str2double(strTmp);
    
    loc0 = strfind(zz.metadataXML, '<CenterY>');
    loc1 = strfind(zz.metadataXML, '</CenterY>');
    strTmp = zz.metadataXML(loc0(1)+9:loc1(1)-1);
    location(1) = str2double(strTmp);
    
    loc0 = strfind(zz.metadataXML, '<Radius>');
    loc1 = strfind(zz.metadataXML, '</Radius>');
    if isempty(loc0)
        % possible ellipse
        loc0 = strfind(zz.metadataXML, '<RadiusX>');
        loc1 = strfind(zz.metadataXML, '</RadiusX>');
        strTmp1 = zz.metadataXML(loc0(1)+8:loc1(1)-1);
        loc0 = strfind(zz.metadataXML, '<RadiusY>');
        loc1 = strfind(zz.metadataXML, '</RadiusY>');
        strTmp2 = zz.metadataXML(loc0(1)+8:loc1(1)-1);
        location(3) = (str2double(strTmp1)+str2double(strTmp2))/2;
    else
        strTmp = zz.metadataXML(loc0(1)+8:loc1(1)-1);
        location(3) = str2double(strTmp);
    end
end

function ablation_region = getAblationRegion(foreground, loc)
    % code from https://www.zhihu.com/question/26551754/answer/33185339
    [x,y,z] = size(foreground);
    num_point = size(loc,1);
    ablation_2d = zeros(x,y);
    ablation_region = zeros(x,y,z);
    for xx = 1:x
        for kk = 1:num_point
            p1 = loc(kk,:);
            if kk < num_point
                p2 = loc(kk+1,:);
            else
                p2 = loc(1,:);
            end
            if p1(1) == p2(1)
                continue;
            elseif xx < min(p1(1),p2(1))
                continue;
            elseif xx > max(p1(1),p2(1))
                continue;
            end
            yy = (xx - p1(1))*(p2(2)-p1(2)) / (p2(1)-p1(1)) + p1(2);
            start_loc = max(ceil(yy),1);
            ablation_2d(xx,start_loc:end) = ablation_2d(xx,start_loc:end) + 1;        
        end
    end
    ablation_2d = mod(ablation_2d,2) == 1;
    for zz = 1:z
        ablation_region(:,:,zz) = ablation_2d.*foreground(:,:,zz);
    end
end
    


function ablation_region = getAblationRegion_Circle(foreground, loc)
    % code from https://www.zhihu.com/question/26551754/answer/33185339
    [x,y,z] = size(foreground);
    ablation_2d = zeros(x,y);
    ablation_region = zeros(x,y,z);
    for xx = 1:x
        for yy = 1:y
            if (xx - loc(1))^2 + (yy - loc(2))^2 <= loc(3)^2
                ablation_2d(xx,yy) = 1;
            end
        end
    end
    for zz = 1:z
        ablation_region(:,:,zz) = ablation_2d.*foreground(:,:,zz);
    end
end
    