function [eig_all, overlay_cl] = principalCv3d(data, sigma)
% use principle curvature to segment connected cells

% get the hessian matrix
[Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, ~] = Hessian3D(data,sigma);

% test each connected component
eig_all = zeros(size(data));
dir_y = zeros(size(data));
dir_x = zeros(size(data));
dir_z = zeros(size(data));

fmap = ones(size(data));
s = regionprops3(fmap, {'VoxelIdxList'});
for i=1:numel(s.VoxelIdxList)
    vox = s.VoxelIdxList{i};
    xx = Dxx(vox); yy = Dyy(vox); zz = Dzz(vox);
    xy = Dxy(vox); xz = Dxz(vox); yz = Dyz(vox);
    
    C = zeros(numel(s.VoxelIdxList{i}),3);
    dir_xyz = zeros(numel(s.VoxelIdxList{i}),3);
    parfor j=1:numel(s.VoxelIdxList{i})
        MM = [xx(j), xy(j), xz(j);...
            xy(j), yy(j), yz(j);...
            xz(j), yz(j), zz(j)];
        [Evec,Eval] = eig(MM);
        dEval = diag(Eval);
        [c,od] = sort(dEval,'descend');
        C(j,:) = c';
        dir_xyz(j,:) = Evec(:, od(1))';
    end
    dir_x(vox) = dir_xyz(:,1);
    dir_y(vox) = dir_xyz(:,2);
    dir_z(vox) = dir_xyz(:,3);
%     eig_all(vox) = C(:,1);
%     eig_all(vox) = abs(C(:,2)) .* (1 - abs(C(:,1))./abs(C(:,2)));
    eig_all(vox) = -C(:,2);
end

end