clc;clear;close all;
dbstop if error

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
addpath ./bfmatlab
xls_name = 'result.xlsx';
xlswrite([file_name '\' xls_name],{'Path'},'Sheet1','A1:A1');
xlswrite([file_name '\' xls_name],{file_name},'Sheet1','B1:B1');

% basic information
czi_file = bfopen([file_name '.czi']);
ome = czi_file{1,4};
rx = double(ome.getPixelsPhysicalSizeX(0).value);
ry = double(ome.getPixelsPhysicalSizeY(0).value);
rz = double(ome.getPixelsPhysicalSizeZ(0).value);
unit_volume = rx*ry*rz;
dv = [rx ry rz];
x = czi_file{1,4}.getPixelsSizeX(0).getValue();
y = czi_file{1,4}.getPixelsSizeY(0).getValue();
z = czi_file{1,4}.getPixelsSizeZ(0).getValue();
c = czi_file{1,4}.getPixelsSizeC(0).getValue();
t = czi_file{1,4}.getPixelsSizeT(0).getValue();
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

xlswrite([file_name '\' xls_name],{'Ablation time'},'Sheet1','A2:A2');
xlswrite([file_name '\' xls_name],time_second(ablation_time),'Sheet1','B2:B2');
xlswrite([file_name '\' xls_name],{'Ablation frame'},'Sheet1','A3:A3');
xlswrite([file_name '\' xls_name],ablation_time,'Sheet1','B3:B3');

% ablation region detection
location = [140.911357 158.143121 139.799631 163.701754 160.088643 168.426593 161.756233 160.088643];
load([file_name '\foreground_c1.mat']);
foreground_c1 = foreground;
[x, y, z] = size(foreground_c1);
a(1) = round(location(2));
a(2) = round(location(1));
b(1) = round(location(4));
b(2) = round(location(3));
c(1) = round(location(6));
c(2) = round(location(5));
d(1) = round(location(8));
d(2) = round(location(7));
c_p = @(x,y) x(1)*y(2) - x(2)*y(1);
ablation_center = (a+b+c+d)/4;
ablation_region = zeros(size(foreground_c1));
for xx = 1:x
    for yy = 1:y
        for zz = 1:z
            p = [xx yy];
            ab = a-b;
            bc = b-c;
            cd = c-d;
            da = d-a;
            ap = a-p;
            bp = b-p;
            cp = c-p;
            dp = d-p;
            if (sign(c_p(ab,ap)) > 0 && sign(c_p(bc,bp)) > 0 && sign(c_p(cd,cp)) > 0 && sign(c_p(da,ap)) > 0) ||...
                    (sign(c_p(ab,ap)) < 0 && sign(c_p(bc,bp)) < 0 && sign(c_p(cd,cp)) < 0 && sign(c_p(da,ap)) < 0)
                if foreground_c1(xx,yy,zz) == 1
                    ablation_region(xx,yy,zz) = 1;
                end
            end
        end
    end
end
[~,~,z_ind] = ind2sub(size(ablation_region),find(ablation_region));
ablation_center = [ablation_center mean(z_ind)];
tic;
distance_c1 = getFrameDist(foreground_c1, dv);
toc

% find 2 min and 10 min
time_minute = (time_second - time_second(ablation_time))/60;
[~, frame_2min] = min(abs(time_minute-2));
[~, frame_6min] = min(abs(time_minute-6));
[~, frame_10min] = min(abs(time_minute-10));
xlswrite([file_name '\' xls_name],{'Frame after 2 min'},'Sheet1','A4:A4');
xlswrite([file_name '\' xls_name],frame_2min,'Sheet1','B4:B4');
xlswrite([file_name '\' xls_name],{'Frame after 6 min'},'Sheet1','A5:A5');
xlswrite([file_name '\' xls_name],frame_6min,'Sheet1','B5:B5');
xlswrite([file_name '\' xls_name],{'Frame after 10 min'},'Sheet1','A6:A6');
xlswrite([file_name '\' xls_name],frame_10min,'Sheet1','B6:B6');
foreground_c1_num = sum(foreground_c1(:));

xlswrite([file_name '\' xls_name],{'Vessel volume'},'Sheet1','A9:A9');
xlswrite([file_name '\' xls_name],{'Glia volume(begining)'},'Sheet1','A10:A10');
xlswrite([file_name '\' xls_name],{'Glia volume(ablation)'},'Sheet1','A11:A11');
xlswrite([file_name '\' xls_name],{'Glia volume(2min)'},'Sheet1','A12:A12');
xlswrite([file_name '\' xls_name],{'Glia volume(6min)'},'Sheet1','A13:A13');
xlswrite([file_name '\' xls_name],{'Glia volume(10min)'},'Sheet1','A14:A14');
xlswrite([file_name '\' xls_name],{'Glia volume(end)'},'Sheet1','A15:A15');
xlswrite([file_name '\' xls_name],foreground_c1_num*unit_volume,'Sheet1','B9:B9');

frame_list = [1 ablation_time frame_2min frame_6min frame_10min t];
for ii = 1:6
    tt = frame_list(ii);
    ind = num2str(1000+tt);
    ind = ind(2:4);
    load([file_name '\foreground_c2\' ind '.mat']);
    foreground_c2 = foreground;
    if ii == 2
        foreground_ablation = foreground;
    end
    foreground_c2_num = sum(foreground_c2(:));
    
    jj = num2str(9+ii);
    xlswrite([file_name '\' xls_name],foreground_c2_num*unit_volume,'Sheet1',['B' jj ':B' jj]);
end
xlswrite([file_name '\' xls_name],{'Total volume of view'},'Sheet1','A16:A16');
xlswrite([file_name '\' xls_name],total_volume,'Sheet1','B16:B16');

%%% soma
ind = num2str(1000+ablation_time);
ind = ind(2:4);
load([file_name '.\soma\' ind '.mat']);
soma_components = bwconncomp(soma);
xlswrite([file_name '\' xls_name],{'Num of cell bodies at ablation'},'Sheet1','A19:A19');
xlswrite([file_name '\' xls_name],soma_components.NumObjects,'Sheet1','B19:B19');
xlswrite([file_name '\' xls_name],{'Cell body ID'},'Sheet1','A20:A20');
xlswrite([file_name '\' xls_name],{'distance2vessel'},'Sheet1','B20:B20');
xlswrite([file_name '\' xls_name],{'distance2ablation_site'},'Sheet1','C20:C20');
for ss = 1:soma_components.NumObjects
    [x_ind, y_ind, z_ind] = ind2sub(size(soma), soma_components.PixelIdxList{ss});
    x_soma = mean(x_ind); y_soma = mean(y_ind); z_soma = mean(z_ind);
    soma2vessel = distance_c1(round(x_soma),round(y_soma),round(z_soma));
    soma2ablation = getP2Vdist([x_soma y_soma z_soma],ablation_center,dv);
    jj = num2str(20+ss);
    xlswrite([file_name '\' xls_name],num2str(ss),'Sheet1',['A' jj ':A' jj]);
    xlswrite([file_name '\' xls_name],soma2vessel,'Sheet1',['B' jj ':B' jj]);
    xlswrite([file_name '\' xls_name],soma2ablation,'Sheet1',['C' jj ':C' jj]);
end
%%% tracking
load([file_name '\data_c2.mat']);
data_all = data2;
load([file_name '\tip_info_reconnect.mat']);
[x, y, z, t] = size(data_all);
addpath ./CINDA

v_max = 12;

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
xlswrite([file_name '\' xls_name],{'track_ID'},'Sheet2','A1:A1');
xlswrite([file_name '\' xls_name],{'first_frame'},'Sheet2','B1:B1');
xlswrite([file_name '\' xls_name],{'last_frame'},'Sheet2','C1:C1');
xlswrite([file_name '\' xls_name],{'total_distance'},'Sheet2','D1:D1');
xlswrite([file_name '\' xls_name],{'end_distance'},'Sheet2','E1:E1');
xlswrite([file_name '\' xls_name],{'tortuosity'},'Sheet2','F1:F1');
xlswrite([file_name '\' xls_name],{'valid_at_ablation'},'Sheet2','G1:G1');
xlswrite([file_name '\' xls_name],{'distance2vessel_at_ablation'},'Sheet2','H1:H1');
xlswrite([file_name '\' xls_name],{'cell_body_identity'},'Sheet2','I1:I1');
%
label_ablation = bwlabeln(foreground_ablation, 26);
for ii = 1:length(trajectories)
    ii
    start_ind = trajectories{ii}(1);
    end_ind = trajectories{ii}(end);
    first_frame = node_table(start_ind,1);
    last_frame = node_table(end_ind,1);
    
    total_distance = 0;
    end_distance = getP2Pdist(node_table(end_ind,2:4), node_table(start_ind,2:4), dv);
    for jj = 2:length(trajectories{ii})
        total_distance = total_distance + getP2Pdist(node_table(trajectories{ii}(jj-1),2:4), node_table(trajectories{ii}(jj),2:4), dv);
    end
    
    t_temp = zeros(length(trajectories{ii}),1);
    for jj = 1:length(trajectories{ii})
        t_temp(jj) = node_table(trajectories{ii}(jj),1);
    end
    if ismember(ablation_time, t_temp)
        valid_at_ablation = 1;
        [~, frame_loc] = ismember(ablation_time, t_temp);
        frame_ind = trajectories{ii}(frame_loc);
        distance2vessel_at_ablation = distance_c1(node_table(frame_ind,2),node_table(frame_ind,3),node_table(frame_ind,4));
        cell_body_identity = 0;
        for jj = 1:soma_components.NumObjects
            if label_ablation(node_table(frame_ind,2),node_table(frame_ind,3),node_table(frame_ind,4)) == label_ablation(soma_components.PixelIdxList{jj}(1))
                cell_body_identity = jj;
            end
        end
    else
        valid_at_ablation = 0;
        distance2vessel_at_ablation = -1;
        cell_body_identity = -1;
    end
 
    
    jj = num2str(ii+1);
    xlswrite([file_name '\' xls_name],ii,'Sheet2',['A' jj ':A' jj]);
    xlswrite([file_name '\' xls_name],first_frame,'Sheet2',['B' jj ':B' jj]);
    xlswrite([file_name '\' xls_name],last_frame,'Sheet2',['C' jj ':C' jj]);
    xlswrite([file_name '\' xls_name],total_distance,'Sheet2',['D' jj ':D' jj]);
    xlswrite([file_name '\' xls_name],end_distance,'Sheet2',['E' jj ':E' jj]);
    xlswrite([file_name '\' xls_name],total_distance/end_distance,'Sheet2',['F' jj ':F' jj]);
    xlswrite([file_name '\' xls_name],valid_at_ablation,'Sheet2',['G' jj ':G' jj]);
    xlswrite([file_name '\' xls_name],distance2vessel_at_ablation,'Sheet2',['H' jj ':H' jj]);
    xlswrite([file_name '\' xls_name],cell_body_identity,'Sheet2',['I' jj ':I' jj]);
end

%% Sheet 3 old variable
xlswrite([file_name '\' xls_name],{'track_ID'},'Sheet3','A1:A1');
xlswrite([file_name '\' xls_name],{'distance_baseline'},'Sheet3','B1:B1');
xlswrite([file_name '\' xls_name],{'speed_baseline'},'Sheet3','C1:C1');
xlswrite([file_name '\' xls_name],{'distance_2min_after_ablation'},'Sheet3','D1:D1');
xlswrite([file_name '\' xls_name],{'speed_2min_after_ablation'},'Sheet3','E1:E1');
xlswrite([file_name '\' xls_name],{'distance_6min_after_ablation'},'Sheet3','F1:F1');
xlswrite([file_name '\' xls_name],{'speed_6min_after_ablation'},'Sheet3','G1:G1');
xlswrite([file_name '\' xls_name],{'distance_10min_after_ablation'},'Sheet3','H1:H1');
xlswrite([file_name '\' xls_name],{'speed_10min_after_ablation'},'Sheet3','I1:I1');
xlswrite([file_name '\' xls_name],{'distance2ablation_start'},'Sheet3','J1:J1');
xlswrite([file_name '\' xls_name],{'distance2ablation_at_ablation'},'Sheet3','K1:K1');
xlswrite([file_name '\' xls_name],{'distance2ablation_end'},'Sheet3','L1:L1');
for ii = 1:length(trajectories)
    ii
    
    t_temp = zeros(length(trajectories{ii}),1);
    for jj = 1:length(trajectories{ii})
        t_temp(jj) = node_table(trajectories{ii}(jj),1);
    end
    start_ind = trajectories{ii}(1);
    end_ind = trajectories{ii}(end);
    if (t_temp(1) < ablation_time) && ismember(ablation_time, t_temp)
        distance_baseline = getTrackDist(trajectories{ii}, node_table, t_temp(1), ablation_time, dv);
        speed_baseline = distance_baseline/(time_second(ablation_time) - time_second(t_temp(1)));
    else
        distance_baseline = -1;
        speed_baseline = -1;
    end
    
    if ismember(ablation_time, t_temp) && ismember(frame_2min, t_temp)
        distance_2min = getTrackDist(trajectories{ii}, node_table, ablation_time, frame_2min, dv);
    else
        distance_2min = -1;
    end
    if ismember(frame_2min, t_temp)
        speed_2min = getSpeed(trajectories{ii}, node_table, frame_2min, dv, time_second);
    else
        speed_2min = -1;
    end
    
    if ismember(ablation_time, t_temp) && ismember(frame_6min, t_temp)
        distance_6min = getTrackDist(trajectories{ii}, node_table, ablation_time, frame_6min, dv);
    else
        distance_6min = -1;
    end
    if ismember(frame_6min, t_temp)
        speed_6min = getSpeed(trajectories{ii}, node_table, frame_6min, dv, time_second);
    else
        speed_6min = -1;
    end
    
    if ismember(ablation_time, t_temp) && ismember(frame_10min, t_temp)
        distance_10min = getTrackDist(trajectories{ii}, node_table, ablation_time, frame_10min, dv);
    else
        distance_10min = -1;
    end
    if ismember(frame_10min, t_temp)
        speed_10min = getSpeed(trajectories{ii}, node_table, frame_10min, dv, time_second);
    else
        speed_10min = -1;
    end
    
    distance2ablation_start = getP2Pdist(node_table(start_ind,2:4), ablation_center, dv);
    distance2ablation_end = getP2Pdist(node_table(end_ind,2:4), ablation_center, dv);    
    if ismember(ablation_time, t_temp)
        [~, frame_loc] = ismember(ablation_time, t_temp);
        frame_ind = trajectories{ii}(frame_loc);
        distance2ablation_ablation = getP2Pdist(node_table(frame_ind,2:4), ablation_center, dv);  
    else
        distance2ablation_ablation = -1;
    end
    
    
    jj = num2str(ii+1);
    xlswrite([file_name '\' xls_name],ii,'Sheet3',['A' jj ':A' jj]);
    xlswrite([file_name '\' xls_name],distance_baseline,'Sheet3',['B' jj ':B' jj]);
    xlswrite([file_name '\' xls_name],speed_baseline,'Sheet3',['C' jj ':C' jj]);
    xlswrite([file_name '\' xls_name],distance_2min,'Sheet3',['D' jj ':D' jj]);
    xlswrite([file_name '\' xls_name],speed_2min,'Sheet3',['E' jj ':E' jj]);
    xlswrite([file_name '\' xls_name],distance_6min,'Sheet3',['F' jj ':F' jj]);
    xlswrite([file_name '\' xls_name],speed_6min,'Sheet3',['G' jj ':G' jj]);
    xlswrite([file_name '\' xls_name],distance_10min,'Sheet3',['H' jj ':H' jj]);
    xlswrite([file_name '\' xls_name],speed_10min,'Sheet3',['I' jj ':I' jj]);
    xlswrite([file_name '\' xls_name],distance2ablation_start,'Sheet3',['J' jj ':J' jj]);
    xlswrite([file_name '\' xls_name],distance2ablation_ablation,'Sheet3',['K' jj ':K' jj]);
    xlswrite([file_name '\' xls_name],distance2ablation_end,'Sheet3',['L' jj ':L' jj]);
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