clc;clear;close all;

file_name = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series';
load([file_name '\data_c2']);
data = data2;
[x,y,z,t] = size(data);
scale = 2;
data_c = zeros(x*scale,y*scale,3,z,t);
% for tt = 1:t
%     for zz = 1:z
%         data_c(:,:,1,zz,tt) = imresize(data(:,:,zz,tt),scale)/2;
%         data_c(:,:,2,zz,tt) = imresize(data(:,:,zz,tt),scale)/2;
%         data_c(:,:,3,zz,tt) = imresize(data(:,:,zz,tt),scale)/2;
%     end
% end
load([file_name '\time_domain_tracking\data_c']);

% sp = [133   141    22]; st = 10;
% sp = [158   254    34]; st = 2;
% sp = [153   209    40]; st = 2;  %bad
% sp = [161   204    40]; st = 8;  %bad
% sp = [174   125    22]; st = 14;
% sp = [191   194     9]; st = 19;
sp = [178   164    35 ]; st = 19;
win = 2;
color = (rand(1,3)/1.5+0.5/1.5)*255;
for tt = st:t
    tt
    load([file_name '\result_gradientDescent\result_gradientDescent_' num2str(tt-1) '_0.01_0.2_1500.mat'])
    vec = [0 0 0];
    for ii = -win:win 
        for jj = -win:win
            for kk = -win/2:win/2
                try
                    vec_temp = [-ux(sp(1)+ii,sp(2)+jj,sp(3)+kk) -uy(sp(1)+ii,sp(2)+jj,sp(3)+kk) -uz(sp(1)+ii,sp(2)+jj,sp(3)+kk)];
                    if norm(vec) < norm(vec_temp)
                        vec = vec_temp;
                        sp_temp = [sp(1)+ii sp(2)+jj sp(3)+kk];
                    end
                catch me
                    disp(me.message);
                end
            end
        end
    end
%     sp = sp_temp;
    vec = round(vec)
    ep = sp_temp + round(vec);
    for ttt = tt:t
        data_c(:,:,:,:,ttt) = draw3Dline(data_c(:,:,:,:,ttt), [sp(1)*scale sp(2)*scale sp(3)], [ep(1)*scale ep(2)*scale ep(3)], 1.5, color);
    end
    sp = ep;
end
save([file_name '\time_domain_tracking\data_c'], 'data_c', '-v7.3');

for tt = 1:t
    ind = num2str(1000+tt);
    ind = ind(2:4);
    tifwrite(uint8(data_c(:,:,:,:,tt)),[file_name '\time_domain_tracking\' num2str(ind)]);
end