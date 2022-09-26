clc;clear;close all;

file_name = '..\SL-092320-slice1-hippo-vessel-Modify Series';
mkdir([file_name '\foreground']);
load([file_name '\data_c2']);

[x,y,z] = size(data);
std = sqrt(max(variance));

x_direction = zeros(1,27);
y_direction = zeros(1,27);
z_direction = zeros(1,27);
ind = 1;
for xx = -1:1
    for yy = -1:1
        for zz = -1:1
            x_direction(ind) = xx;
            y_direction(ind) = yy;
            z_direction(ind) = zz;
            ind = ind + 1;
        end
    end
end

% 10 neighbor connecivty
conn = cat(3, [0 0 0; 0 1 0; 0 0 0],...
              [1 1 1; 1 1 1; 1 1 1],...
              [0 0 0; 0 1 0; 0 0 0]);
        
% for tt = 1:1
%     tt 
% ind = num2str(1000+tt); 
% data = load(['.\StabilizedData\' ind(2:4) '.mat']);
% data = data.temp;
    
components_cell = cell(25,2);

tic;
for thre = 1:25
    thre
    binaryData = zeros(x,y,z);
    binaryData(data<thre*10) = 0;
    binaryData(data>=thre*10) = 1;
    
    components = bwconncomp(binaryData,conn);
    components_len = length(components.PixelIdxList);
    z_scores = zeros(components_len,1);
    for com_ind = 1:components_len
%        com_ind       
       com_element = components.PixelIdxList{com_ind};
       if length(com_element)>= 30
%            time1=toc;
%            element_graph = zeros(x,y,z);
%            neighbor_graph = zeros(x,y,z);
%            [x_location,y_location,z_location] = ind2sub([x,y,z],com_element);
%            
%            for element_ind = 1:length(com_element)
%                element_graph(x_location(element_ind),y_location(element_ind),z_location(element_ind)) = 1;
%                neighbor_graph(max(1,x_location(element_ind)-1):min(x_location(element_ind)+1,x),...
%                    max(1,y_location(element_ind)-1):min(y_location(element_ind)+1,y),max(1,z_location(element_ind)-1):min(z_location(element_ind)+1,z)) = 1;
%            end 

            
           [x_position,y_position,z_position] = ind2sub([x,y,z],com_element);
           com_neighbor = [];
           for dir_ind = 1:27
              x_neighbor = min(max(x_position + x_direction(dir_ind),1),x); 
              y_neighbor = min(max(y_position + y_direction(dir_ind),1),y); 
              z_neighbor = min(max(z_position + z_direction(dir_ind),1),z); 
              pix1 = sub2ind([x,y,z],x_neighbor,y_neighbor,z_neighbor);
              pix1 = pix1(~ismember(pix1,com_element));
              pix1 = pix1(~ismember(pix1,com_neighbor));
              com_neighbor = [com_neighbor;pix1];
           end
           M = numel(com_element);
           N = numel(com_neighbor);
           L = mean(data(com_element)) - mean(data(com_neighbor));

%            neighbor_graph = neighbor_graph&(1-element_graph);
%            M = sum(element_graph(:));
%            N = sum(neighbor_graph(:));
%            element_data = element_graph.*data;
%            neighbor_data = neighbor_graph.*data;
%            L = sum(element_data(:))/M - sum(neighbor_data(:))/N;  
           pdfz = normpdf(norminv(N/(M+N),0,1));
           cdfz = normcdf(norminv(N/(M+N),0,1));
           mu = pdfz/(1-cdfz) + pdfz/cdfz;
           sigma = sqrt((1+norminv(N/(M+N),0,1)*pdfz/(1-cdfz)-(pdfz/(1-cdfz))^2)/M + (1-norminv(N/(M+N),0,1)*pdfz/cdfz-(pdfz/cdfz)^2)/N);
           z_scores(com_ind) = (L/std - mu)/sigma;

%            if M<=1000 && N<=1000
%                z_scores(com_ind) = (L - Mu(M,N))/Sigma(M,N);
%            else
%                z_scores(com_ind) = (L - integralMu(M,N,1))/integralSigma(M,N,1);
%            end  
%          z_scores(com_ind)
%            if z_scores(com_ind) > 2000
%                foreground = 1- (1-foreground).*(1-element_graph);
%            end
       end
    end
    components_cell(thre,1) = {components.PixelIdxList};
    components_cell(thre,2) = {z_scores};
end
toc
save([file_name '\foreground'],'components_cell');
    