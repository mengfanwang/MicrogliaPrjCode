clc;clear;close all;

file_name = '091919 slice1 cortex 4mM isoflurane';

mkdir([file_name '\result_post']);
load([file_name '\matchingResult']);
load([file_name '\graph']);
load([file_name '\result']);
load([file_name '\cell']);
load([file_name '\cell_bg']);
result_all = result;

cell_node = median(cellfun(@length,Node));
cell_volume = zeros(1,length(Node));
for tt = 1:length(Node)
    cell_volume(tt) = sum(arrayfun(@(x)length(x.voxel),Node{tt}));
end
cell_volume = median(cell_volume);

tic;
% find unmapped components
for tt = 1:length(Node)-1

[g,pos1,pos2,s,t,s_sep,t_sep,s_com{tt},s_long] = buildGraph(Node{tt},Node{tt+1},G{tt},G{tt+1},matchingResult{tt});
com_label = ones(1,length(s_com{tt}));
for ii = 1:length(s_com{tt})
    if length(s_com{tt}{ii}) >=4
        com_label(ii) = 0;
    end
end
s_com{tt}(logical(com_label)) = [];

% figure(1); 
% h1 = plot(g,'XData',[pos1(1,:) pos2(1,:)],'YData',[pos1(2,:) pos2(2,:)],'ZData',[pos1(3,:) pos2(3,:)]);%,'NodeLabel',nLabel);%,'Layout','force');
% highlight(h1,s,t,'EdgeColor','g');
% highlight(h1,[s_sep t_sep],'NodeColor','g');
% highlight(h1,s_long,'NodeColor','r');
% set(gca,'Color','k');

% single
% marker_size = zeros(1,length(Node{tt}));
% for ii = 1:length(Node{tt})
%     marker_size(ii) = sqrt(length(Node{tt}(ii).voxel))/5;
% end
% figure(1); 
% h1 = plot(G{tt},'XData',pos1(1,:),'YData',pos1(2,:),'ZData',pos1(3,:));%,'NodeLabel',nLabel);%,'Layout','force');
% h1.MarkerSize = marker_size;
% set(gca,'Color','k');

[g,pos1,pos2,s,t,s_sep,t_sep,t_com{tt+1},t_long] = buildGraph(Node{tt+1},Node{tt},G{tt+1},G{tt},matchingResult{tt}');
com_label = ones(1,length(t_com{tt+1}));
for ii = 1:length(t_com{tt+1})
    if length(t_com{tt+1}{ii}) >=4
        com_label(ii) = 0;
    end
end
t_com{tt+1}(logical(com_label)) = [];

% figure(2); 
% h1 = plot(g,'XData',[pos1(1,:) pos2(1,:)],'YData',[pos1(2,:) pos2(2,:)],'ZData',[pos1(3,:) pos2(3,:)]);%,'NodeLabel',nLabel);%,'Layout','force');
% highlight(h1,s,t,'EdgeColor','g');
% % highlight(h1,[s_sep t_sep],'NodeColor','g');
% highlight(h1,t_long,'NodeColor','r');
% set(gca,'Color','k');
end

% remove big components
for tt = 1:length(Node)
    ind = num2str(1000+tt);
    ind = ind(2:4);
    result = result_all{tt};
    [x_size, y_size, ~, z_size] = size(result);
    over_connection{tt} = zeros(x_size,y_size,z_size);
    if tt < 52
        for ii = 1:length(s_com{tt})
            if length(s_com{tt}{ii}) > cell_node/2
                for jj = 1:length(s_com{tt}{ii})
                    over_connection{tt}(Node{tt}(s_com{tt}{ii}(jj)).voxel) = 1;
                end
            end
        end
    end
    if tt > 1
        for ii = 1:length(t_com{tt})
            if length(t_com{tt}{ii}) > cell_node/2
                for jj = 1:length(t_com{tt}{ii})
                    over_connection{tt}(Node{tt}(t_com{tt}{ii}(jj)).voxel) = 1;
                end
            end
        end    
    end
%     for zz = 1:z_size
%         result(:,:,1,zz) = over_connection(:,:,zz).*128 + result(:,:,1,zz).*(1-over_connection(:,:,zz));
%         result(:,:,2,zz) = over_connection(:,:,zz).*128 + result(:,:,2,zz).*(1-over_connection(:,:,zz));
%         result(:,:,3,zz) = over_connection(:,:,zz).*128 + result(:,:,3,zz).*(1-over_connection(:,:,zz));
%     end
    
    cell{tt} = cell{tt} - over_connection{tt};
end
% find small components
for tt = 1:length(Node)
    tt
    ind = num2str(1000+tt);
    ind = ind(2:4);
    
    fore_reg = cell_bg{tt} - cell{tt};
    
    fore_com = bwconncomp(fore_reg,26);
    fore_com.PixelIdxList(cellfun(@length,fore_com.PixelIdxList) > cell_volume/2) = [];
    fore_new = zeros(size(fore_reg));
    
    if tt > 1
        s_com{tt-1}(cellfun(@length,s_com{tt-1}) > cell_node/2) = [];
        for ii = 1:length(s_com{tt-1})
            for jj = 1:length(s_com{tt-1}{ii})
                for kk = 1:length(fore_com.PixelIdxList)
                    if any(ismember(Node{tt-1}(s_com{tt-1}{ii}(jj)).voxel,fore_com.PixelIdxList{kk}))
                        fore_new(fore_com.PixelIdxList{kk}) = 1;
                    end
                end
            end
        end
    end
    if tt < 52
        t_com{tt+1}(cellfun(@length,t_com{tt+1}) > cell_node/2) = [];
        for ii = 1:length(t_com{tt+1})
            for jj = 1:length(t_com{tt+1}{ii})
                for kk = 1:length(fore_com.PixelIdxList)
                    if any(ismember(Node{tt+1}(t_com{tt+1}{ii}(jj)).voxel,fore_com.PixelIdxList{kk}))
                        fore_new(fore_com.PixelIdxList{kk}) = 1;
                    end
                end
            end
        end
    end

    result = result_all{tt};
    [x_size, y_size, ~, z_size] = size(result);
    for zz = 1:z_size
        result(:,:,1,zz) = over_connection{tt}(:,:,zz).*128 + result(:,:,1,zz).*(1-over_connection{tt}(:,:,zz) - fore_new(:,:,zz));
        result(:,:,2,zz) = over_connection{tt}(:,:,zz).*128 + result(:,:,2,zz).*(1-over_connection{tt}(:,:,zz) - fore_new(:,:,zz)) + fore_new(:,:,zz)*255;
        result(:,:,3,zz) = over_connection{tt}(:,:,zz).*128 + result(:,:,3,zz).*(1-over_connection{tt}(:,:,zz) - fore_new(:,:,zz));
    end
    tifwrite(uint8(result),[file_name '\result_post\' ind]);
end


function [g,pos1,pos2,s,t,s_sep,t_sep,s_com,s_long] = buildGraph(Node1,Node2,G1,G2,map)
    n1 = length(Node1); n2 = length(Node2);
    m1 = height(G1.Edges); m2 = height(G2.Edges);

    pos1 = zeros(3,n1); pos2 = zeros(3,n2);
    for ii = 1:n1
        pos1(1,ii) = Node1(ii).x;
        pos1(2,ii) = Node1(ii).y;
        pos1(3,ii) = Node1(ii).z;
    end
    for ii = 1:n2
        pos2(1,ii) = Node2(ii).x;
        pos2(2,ii) = Node2(ii).y;
        pos2(3,ii) = Node2(ii).z;
    end
    Edge = [G1.Edges.EndNodes; G2.Edges.EndNodes + n1];
    g = digraph(Edge(:,1),Edge(:,2));

    [s,t] = find(map);
    dis = zeros(1,length(s));
    for ii = 1:length(s)
        v = [Node2(t(ii)).x - Node1(s(ii)).x  Node2(t(ii)).y - Node1(s(ii)).y  Node2(t(ii)).z - Node1(s(ii)).z];
        dis(ii) = norm(v);
    end
    s = s(dis<5);
    s_sep = 1:n1;
    s_sep(s) = [];
    
    t = t(dis<5);
    t_sep = 1:n2;
    t_sep(t) = [];
    t_sep = t_sep + n1;

    s_com = findSepCom(n1,s_sep,Node1);

    s_long = [];
    for ii = 1:length(s_com)
        if length(s_com{ii}) >= 4
            s_long = [s_long s_com{ii}];
        end
    end

    t = t + n1;
    g = addedge(g,s,t);
end

function s_com = findSepCom(n1,s_sep,Frame)
    s_label = zeros(1,n1);
    s_label(s_sep) = 1;
    com_ind = 0;
    s_com = {};
    for ii = 1:length(s_sep)
        if s_label(s_sep(ii)) == 1
            com_ind = com_ind + 1;
            s_com{com_ind} = [];
            candidate_list = s_sep(ii);
            while candidate_list
                candidate = candidate_list(1);
                candidate_list(1) = [];
                if s_label(candidate) == 1
                    s_label(candidate) = 0;
                    s_com{com_ind} = [s_com{com_ind} candidate];
                    candidate_list = [candidate_list Frame(candidate).parent Frame(candidate).child];
                end
            end
        end
    end
end