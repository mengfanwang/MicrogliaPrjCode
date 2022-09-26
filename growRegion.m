function [pix0,candidate_pix,neighbor_record,summation] = growRegion(pix0,dF,candidate_pix,filter,neighbor_record,summation,sigma,minIntesity)

[H,W,T] = size(dF);

if(~exist('minIntesity','var'))
    minIntesity = 2*sigma;
end

N = sum(filter(:));

n_neighbors = double(neighbor_record(candidate_pix));
predicted = summation(candidate_pix);
predicted = predicted./n_neighbors;

%% judge
prior1 =n_neighbors/N;
prior0 = 1 - prior1;
thresholds = min(predicted,max(max(minIntesity,predicted/10),predicted/2 - sigma^2./predicted.*log(prior1./prior0)));
selected = dF(candidate_pix)>thresholds;
add_pixels = candidate_pix(selected);

%% update
pix0 = [pix0;add_pixels];
dw = zeros(27,1);
dh = zeros(27,1);
dz = zeros(27,1);
cnt = 1;
for x = -1:1
    for y = -1:1
        for z = -1:1
            dw(cnt) = x;
            dh(cnt) = y;
            dz(cnt) = z;
            cnt = cnt + 1;
        end
    end
end

candidate_pix = [];
for kk = 1:27
    [ih0,iw0,it0] = ind2sub([H,W,T],add_pixels);
    ih1 = max(min(ih0+dh(kk),H),1);
    iw1 = max(min(iw0+dw(kk),W),1);
    it1 = max(min(it0+dz(kk),T),1);
    pix1 = sub2ind([H,W,T],ih1,iw1,it1);
    candidate_pix = [pix1;candidate_pix];
end
candidate_pix = setdiff(candidate_pix,pix0);
candidate_pix = unique(candidate_pix);

%% update matrix
[dx,dy,dz] = ind2sub(size(filter),find(filter>0));
radius = (size(filter,1)-1)/2;
dx = dx - radius - 1;
dy = dy - radius - 1;
dz = dz - radius - 1;

intensities = dF(add_pixels);
for i = 1:numel(dx)
    [ih,iw,it] = ind2sub([H,W,T],add_pixels);
    ih0 = min(max(ih + dx(i),1),H);
    iw0 = min(max(iw + dy(i),1),W);
    it0 = min(max(it + dz(i),1),T);
    add_pix_shift = sub2ind([H,W,T],ih0,iw0,it0);
    summation(add_pix_shift) = summation(add_pix_shift) + intensities;
    neighbor_record(add_pix_shift) = neighbor_record(add_pix_shift)+1;
end




end