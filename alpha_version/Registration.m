clc;clear;close all;

reg_label = 1;
file_name = '091919 slice1 cortex 4mM isoflurane';

mkdir([file_name '\foreground_reg']);
mkdir([file_name '\soma_reg']);
load([file_name '\foreground']);
load([file_name '\soma']);
foreground_all = foreground;
soma_all = soma;

diff = zeros(size(foreground_all,4)-1,3);
if reg_label
    for tt = 1:size(foreground_all,4)-1

        soma1 = soma(:,:,:,tt);
        soma2 = soma(:,:,:,tt+1);

        component1 = bwconncomp(soma1,26);
        component1_len = length(component1.PixelIdxList);
        component2 = bwconncomp(soma2,26);
        component2_len = length(component2.PixelIdxList);

        info = [];  % 8 column: x1 x2 y1 y2 z1 z2 com1 com2

        for ii = 1:component1_len
            for jj = 1:component2_len
                if length(component1.PixelIdxList{ii}) > 100 && length(component2.PixelIdxList{jj}) > 100 % remove small regions
                if any(ismember(component1.PixelIdxList{ii},component2.PixelIdxList{jj}))
                    [x_ind1,y_ind1,z_ind1] = ind2sub(size(soma1),component1.PixelIdxList{ii});
                    [x_ind2,y_ind2,z_ind2] = ind2sub(size(soma2),component2.PixelIdxList{jj});
                    info = [info; mean(x_ind1) mean(x_ind2) mean(y_ind1) mean(y_ind2) mean(z_ind1) mean(z_ind2) ii jj];
                end
                end
            end
        end

        diff(tt,:) = [mean(info(:,2)-info(:,1)) mean(info(:,4)-info(:,3)) mean(info(:,6)-info(:,5))];
    end
end

diff = -cumsum(diff);
diff = [0 0 0;diff];
left_boundary = [floor(min(diff(:,1))) floor(min(diff(:,2))) floor(min(diff(:,3)))];
right_boundary = [ceil(max(diff(:,1))) ceil(max(diff(:,2))) ceil(max(diff(:,3)))];
[x_size,y_size,z_size] = size(soma(:,:,:,1));
foreground_reg = zeros(x_size+right_boundary(1)-left_boundary(1), y_size+right_boundary(2)-left_boundary(2),...
                    z_size+right_boundary(3)-left_boundary(3),size(foreground_all,4));
soma_reg = foreground_reg;
for tt = 1:size(foreground_all,4)
    ind = num2str(1000+tt);
    ind = ind(2:4);
    
    foreground = foreground_all(:,:,:,tt);
    foreground_temp = zeros(x_size+right_boundary(1)-left_boundary(1), y_size+right_boundary(2)-left_boundary(2),...
                        z_size+right_boundary(3)-left_boundary(3));
    foreground_temp(round(diff(tt,1))-left_boundary(1)+1:round(diff(tt,1))-left_boundary(1)+x_size,...
             round(diff(tt,2))-left_boundary(2)+1:round(diff(tt,2))-left_boundary(2)+y_size,...
             round(diff(tt,3))-left_boundary(3)+1:round(diff(tt,3))-left_boundary(3)+z_size) = foreground;
    foreground_reg(:,:,:,tt) = foreground_temp;
    
    soma = soma_all(:,:,:,tt);
    soma_temp = zeros(x_size+right_boundary(1)-left_boundary(1), y_size+right_boundary(2)-left_boundary(2),...
                        z_size+right_boundary(3)-left_boundary(3));
    soma_temp(round(diff(tt,1))-left_boundary(1)+1:round(diff(tt,1))-left_boundary(1)+x_size,...
             round(diff(tt,2))-left_boundary(2)+1:round(diff(tt,2))-left_boundary(2)+y_size,...
             round(diff(tt,3))-left_boundary(3)+1:round(diff(tt,3))-left_boundary(3)+z_size) = soma;
    soma_reg(:,:,:,tt) = soma_temp;
    
    
    tifwrite(uint8(foreground_temp*255),[file_name '\foreground_reg\' ind]);
    tifwrite(uint8(soma_temp*255),[file_name '\soma_reg\' ind]);
end
foreground_reg = logical(foreground_reg);
save([file_name '\foreground_reg'],'foreground_reg');
soma_reg = logical(soma_reg);
save([file_name '\soma_reg'],'soma_reg');
