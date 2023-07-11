function [vidOut, vIdx, edge_Flag, loc_org_xyz] = crop3D(vidIn, yxz, shift)
% crop a region from vidIn
%INPUT:
% vidIn: the vid to crop from
% yxz: different form. 1, n-by-1 or n-by-3 matix, which indicates the 
% overall voxel index; 2, TODO
% shift: (l+2*shift(1))*(w+2*shift(2))*(h+2*shift(shift(1)))

%OUTPUT:
% vidOut: the cropped region

% contact: ccwang@vt.edu
edge_Flag = false;
if size(yxz,2) ==1
    [yy,xx,zz] = ind2sub(size(vidIn),yxz);
end
if length(shift)==1
    shift = [shift shift shift];
end
[h,w,z] = size(vidIn);

ymin = max(1, min(yy)-shift(1)); ymax = min(h, max(yy)+shift(1));
xmin = max(1, min(xx)-shift(2)); xmax = min(w, max(xx)+shift(2));
zmin = max(1, min(zz)-shift(3)); zmax = min(z, max(zz)+shift(3));

loc_org_xyz = [xmin, ymin, zmin];
if min(yy) == 1 || max(yy) == h || min(xx) == 1 || max(xx) == w% || min(zz) == 1 || max(zz) == z
    edge_Flag = true;
end
vidOut = vidIn(ymin:ymax, xmin:xmax, zmin:zmax);

vIdx = [];
if nargout > 1
    [R, C, P] = ndgrid(ymin:ymax, xmin:xmax, zmin:zmax);
    vIdx = sub2ind(size(vidIn), R, C, P);
end
end