clc;clear;close all;
dbstop if error

%%
% dataset = {'HA1','HA2','HA3','HA4','SL1','SL2','SL3','SL4'};
% for dd = 1:8
path = ['D:\dropbox\Fritz\6 files for In-Vivo analysis\'];
files = dir(path);

for ii = 3:length(files)
    ii
    if strcmp(files(ii).name(end-3:end), '.czi')
        file_name = [path files(ii).name(1:end-4)];
        file_name
        ReadFile(file_name);
%         Registration(file_name);
%         VarianceStabilization(file_name);
%         vessel_OrderStatisticThresholding_v6(file_name);
%         glia_OrderStatisticThresholding_v6(file_name);
%         SomaDetection(file_name);
%         distance2probability(file_name);
%         distance2Maximum(file_name);  
%         Output(file_name, path);
%         Output_microglia_ablation(file_name);
    
    end
end

% end

%%
dbstop if error
path = 'D:\dropbox\Modify Series Data-3nd\Priority 4 files\';
file_name = [path 'SL-092120-slice1-cortex-microglia-Modify Series'];
file_name
Output(file_name, path);