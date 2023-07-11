clc;clear;close all;
dbstop if error

file_name = 'D:\dropbox\Modify Series Data-3nd\Priority 1 files\060221 slice1 hippocampus before vessel_Modify Series';
addpath ./bfmatlab
% get location
zz = czifinfo([file_name '.czi']); 
loc1 = strfind(zz.metadataXML, '<Points>');
loc2 = strfind(zz.metadataXML, '</Points>');
strTmp = zz.metadataXML(loc1(1)+8:loc2(1)-1);
loc1 = strfind(strTmp, ' ');
loc2 = strfind(strTmp, ',');
loc = sort([loc1 loc2]);
location = zeros(1, 8);
location(1) = str2num(strTmp(1:loc(1)-1));
for ii = 2:7
    location(ii) = str2num(strTmp(loc(ii-1)+1:loc(ii)-1));
end
location(8) = str2num(strTmp(loc(7)+1:end));

load([file_name '\registration_c1.mat']);
x_min = floor(min(arrayfun(@(x) x.shifts(1), shifts1)));
y_min = floor(min(arrayfun(@(x) x.shifts(2), shifts1)));
x_bias = shifts1(1).shifts(1) - x_min;
y_bias = shifts1(1).shifts(2) - y_min;
location = location + [y_bias x_bias y_bias x_bias y_bias x_bias y_bias x_bias];

% plot abalation site in 3D
load([file_name '\foreground_c1.mat']);
[x, y, z] = size(foreground);
a(1) = round(location(2));
a(2) = round(location(1));
b(1) = round(location(4));
b(2) = round(location(3));
c(1) = round(location(6));
c(2) = round(location(5));
d(1) = round(location(8));
d(2) = round(location(7));
c_p = @(x,y) x(1)*y(2) - x(2)*y(1);
fore2 = zeros(size(foreground));
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
                if foreground(xx,yy,zz) == 1
                    fore2(xx,yy,zz) = 1;
                end
            end
        end
    end
end
im = zeros(x,y,3,z);
for zz = 1:z
    im(:,:,1,zz) = foreground(:,:,zz)+fore2(:,:,zz);
    im(:,:,2,zz) = foreground(:,:,zz) - fore2(:,:,zz);
    im(:,:,3,zz) = foreground(:,:,zz) - fore2(:,:,zz);
end
tifwrite(uint8(im*255),[file_name '\foreground_c1']);




