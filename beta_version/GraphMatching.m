clc;clear;close all;

file_name = '091919 slice1 cortex 4mM isoflurane';

addpath CINDA
load([file_name '\graph']);
v_max = 5;

num_frame = cellfun(@length,Node);
detection_arcs = zeros(sum(num_frame),4);
node_xyz = zeros(sum(num_frame),3);
node_table = zeros(sum(num_frame),3);
transition_arcs = zeros(sum(num_frame(1:end-1).*num_frame(2:end)),3);
transition_ind = 1;
for ii = 1:length(num_frame)
    ii
    for jj = 1:num_frame(ii)
        ind = sum(num_frame(1:ii-1))+jj;
        detection_arcs(ind,:) = [ind v_max/2 v_max/2 -v_max];
        node_xyz(ind,:) = [Node{ii}(jj).x Node{ii}(jj).y Node{ii}(jj).z];
        node_table(ind,:) = [ind ii jj];
        if ii~= length(num_frame)
            for kk = 1:num_frame(ii+1)
                transition_arcs(transition_ind,:) = [ind sum(num_frame(1:ii))+kk ...
                    sqrt((Node{ii}(jj).x-Node{ii+1}(kk).x)^2 + (Node{ii}(jj).y-Node{ii+1}(kk).y)^2 + (Node{ii}(jj).z-Node{ii+1}(kk).z)^2)];
                transition_ind = transition_ind + 1;
            end
        end
    end
end

[trajectories, costs] = mcc4mot(detection_arcs,transition_arcs);
trajectories(cellfun(@length,trajectories)==1) = [];

% build matching result
for tt = 1:51
    matchingResult{tt} = zeros(num_frame(tt),num_frame(tt+1));
end
for ii = 1:length(trajectories)
    node_table_temp = node_table(trajectories{ii},:);
    for jj = 1:length(trajectories{ii}) - 1
        matchingResult{node_table_temp(jj,2)}(node_table_temp(jj,3),node_table_temp(jj+1,3)) = 1;
    end
end

save([file_name '\matchingResult'],'matchingResult');