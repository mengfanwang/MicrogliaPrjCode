% clc;clear;close all;
% s.vrtx = [1 1 0; 2 1 0; 1 2 0; 2 2 0;];
% s.nvrtx = 4;
% % s = S3D(s)

% % for ii = 1:10000
% s.vrtx = [1 1 0; 2 1 0; 1 2 0; 2 2 0];
% % s.vrtx = rand(4,3);
% s.vrtx
% s.nvrtx = 4;
% s1 = s;
% save('s1');
% % load s1;
% % s = s1;
% dlmwrite('C:\Users\mengf\Desktop\openGJK_v2\userP.dat',4);
% dlmwrite('C:\Users\mengf\Desktop\openGJK_v2\userP.dat',s.vrtx,'delimiter',' ','precision',15,'-append');
% s = S3D(s);
% 
% if s.nvrtx == 1
%     distance = norm(s.lambda(1)*s.vrtx(1,:));
% else
%     distance = norm(sum(s.lambda(1:s.nvrtx)'.*s.vrtx(1:s.nvrtx,:)));
% end
% [x,distance2] = findDistance(s.vrtx');
% if abs(distance- distance2)> 0.01
%     error('da')
% end
% distance
% % end

function s = subalgorithm(s)

% s.vrtx = line_ind';
% s.nvrtx = size(line_ind,2);
if s.nvrtx == 4
    s = S3D(s);
elseif s.nvrtx == 3
    s = S2D(s);
elseif s.nvrtx == 2
    s = S1D(s);
elseif s.nvrtx == 1
    s.lambda = 1;
end

% Wk = s.vrtx(1:s.nvrtx,:)';
% lambda = s.lambda;
% if s.nvrtx == 1
%     distance = norm(s.lambda(1)*s.vrtx(1,:));
% else
%     distance = norm(sum(s.lambda(1:s.nvrtx)'.*s.vrtx(1:s.nvrtx,:)));
% end

end

function s = S3D(s)
% S3D version 2
% vrtx = s.vrtx(1:4,:)';
% for ii = 1:4
%     m_temp = vrtx;
%     m_temp(:,ii) = [];
%     C(ii) = (-1)^ii*det(m_temp);
% end

C(4) =  s.vrtx(1,1)*s.vrtx(2,2)*s.vrtx(3,3) + s.vrtx(2,1)*s.vrtx(3,2)*s.vrtx(1,3) + s.vrtx(3,1)*s.vrtx(1,2)*s.vrtx(2,3)...
        -s.vrtx(1,3)*s.vrtx(2,2)*s.vrtx(3,1) - s.vrtx(2,3)*s.vrtx(3,2)*s.vrtx(1,1) - s.vrtx(3,3)*s.vrtx(1,2)*s.vrtx(2,1);
C(3) = -s.vrtx(1,1)*s.vrtx(2,2)*s.vrtx(4,3) - s.vrtx(2,1)*s.vrtx(4,2)*s.vrtx(1,3) - s.vrtx(4,1)*s.vrtx(1,2)*s.vrtx(2,3)...
        +s.vrtx(1,3)*s.vrtx(2,2)*s.vrtx(4,1) + s.vrtx(2,3)*s.vrtx(4,2)*s.vrtx(1,1) + s.vrtx(4,3)*s.vrtx(1,2)*s.vrtx(2,1);
C(2) =  s.vrtx(1,1)*s.vrtx(3,2)*s.vrtx(4,3) + s.vrtx(3,1)*s.vrtx(4,2)*s.vrtx(1,3) + s.vrtx(4,1)*s.vrtx(1,2)*s.vrtx(3,3)...
        -s.vrtx(1,3)*s.vrtx(3,2)*s.vrtx(4,1) - s.vrtx(3,3)*s.vrtx(4,2)*s.vrtx(1,1) - s.vrtx(4,3)*s.vrtx(1,2)*s.vrtx(3,1);
C(1) = -s.vrtx(2,1)*s.vrtx(3,2)*s.vrtx(4,3) - s.vrtx(3,1)*s.vrtx(4,2)*s.vrtx(2,3) - s.vrtx(4,1)*s.vrtx(2,2)*s.vrtx(3,3)...
        +s.vrtx(2,3)*s.vrtx(3,2)*s.vrtx(4,1) + s.vrtx(3,3)*s.vrtx(4,2)*s.vrtx(2,1) + s.vrtx(4,3)*s.vrtx(2,2)*s.vrtx(3,1);
  
detM = sum(C);

FacesTest(1) = CompareSigns(detM,C(1));
FacesTest(2) = CompareSigns(detM,C(2));
FacesTest(3) = CompareSigns(detM,C(3));
FacesTest(4) = CompareSigns(detM,C(4));


if all(FacesTest)
    s.lambda = C/detM;
else
    d_ = inf;
    for jj = 1:4
%         if SAMESIGN(detM,-C(jj))
        if ~FacesTest(jj)
            stmp = s;
            stmp.vrtx(jj,:) = [];
            stmp = S2D(stmp);
            if stmp.nvrtx == 1
                d = norm(stmp.lambda(1)*stmp.vrtx(1,:));
            else
                d = norm(sum(stmp.lambda(1:stmp.nvrtx)'.*stmp.vrtx(1:stmp.nvrtx,:)));
            end
        if d < d_
            smin = stmp;
            d_ = d;
        end    
        end
    end
    s = smin;
end
end





function s = S2D(s)
% S2D version 2
% vrtx = s.vrtx(1:3,:)';
% s1 = vrtx(:,1);
% s2 = vrtx(:,2);
% s3 = vrtx(:,3);
u = s.vrtx(2,:) - s.vrtx(1,:);
v = s.vrtx(3,:) - s.vrtx(1,:);
n(1) = u(2)*v(3) - u(3)*v(2);
n(2) = u(3)*v(1) - u(1)*v(3);
n(3) = u(1)*v(2) - u(2)*v(1);
% n = cross(s2-s1,s3-s1);

p = (s.vrtx(1,:)*n')*n'/(n*n');

% for ii = 1:3
%     m_temp = vrtx;
%     m_temp = [m_temp; 1 1 1;];
%     m_temp(ii,:) = [];
%     nu_test(ii) = det(m_temp);
% end

nu_test(1) = s.vrtx(1,2)*s.vrtx(2,3) + s.vrtx(2,2)*s.vrtx(3,3) + s.vrtx(3,2)*s.vrtx(1,3) ...
            - s.vrtx(1,3)*s.vrtx(2,2) - s.vrtx(2,3)*s.vrtx(3,2) - s.vrtx(3,3)*s.vrtx(1,2);
nu_test(2) = s.vrtx(1,1)*s.vrtx(2,3) + s.vrtx(2,1)*s.vrtx(3,3) + s.vrtx(3,1)*s.vrtx(1,3) ...
            - s.vrtx(1,3)*s.vrtx(2,1) - s.vrtx(2,3)*s.vrtx(3,1) - s.vrtx(3,3)*s.vrtx(1,1);
nu_test(3) = s.vrtx(1,1)*s.vrtx(2,2) + s.vrtx(2,1)*s.vrtx(3,2) + s.vrtx(3,1)*s.vrtx(1,2) ...
            - s.vrtx(1,2)*s.vrtx(2,1) - s.vrtx(2,2)*s.vrtx(3,1) - s.vrtx(3,2)*s.vrtx(1,1);
        
nu_fabs = abs(nu_test);
[~,indexJ] = max(nu_fabs);
umax = nu_test(indexJ);

indexI = [1 2 3];
indexI(indexJ) = [];
C(1) = p(indexI(1))*s.vrtx(2,indexI(2)) + s.vrtx(2,indexI(1))*s.vrtx(3,indexI(2)) + s.vrtx(3,indexI(1))*p(indexI(2)) ...
     - p(indexI(2))*s.vrtx(2,indexI(1)) - s.vrtx(2,indexI(2))*s.vrtx(3,indexI(1)) - s.vrtx(3,indexI(2))*p(indexI(1));
C(2) = s.vrtx(1,indexI(1))*p(indexI(2)) + p(indexI(1))*s.vrtx(3,indexI(2)) + s.vrtx(3,indexI(1))*s.vrtx(1,indexI(2)) ...
     - s.vrtx(1,indexI(2))*p(indexI(1)) - p(indexI(2))*s.vrtx(3,indexI(1)) - s.vrtx(3,indexI(2))*s.vrtx(1,indexI(1));
C(3) = s.vrtx(1,indexI(1))*s.vrtx(2,indexI(2)) + s.vrtx(2,indexI(1))*p(indexI(2)) + p(indexI(1))*s.vrtx(1,indexI(2)) ...
     - s.vrtx(1,indexI(2))*s.vrtx(2,indexI(1)) - s.vrtx(2,indexI(2))*p(indexI(1)) - p(indexI(2))*s.vrtx(1,indexI(1));




% vrtx(indexJ,:) = [];
% vrtx = [vrtx; 1 1 1;];
% p(indexJ) = [];
% p = [p; 1];
% indexI = [1 2 3];
% indexI(indexJ) = [];
% B = vrtx\p*umax;

FacesTest(1) = CompareSigns(umax,C(1));
FacesTest(2) = CompareSigns(umax,C(2));
FacesTest(3) = CompareSigns(umax,C(3));
if all(FacesTest)
    s.lambda = C/umax;
    s.nvrtx = 3;
else
    d_ = inf;
    for jj = 1:3
%         if SAMESIGN(umax,-B(jj))
        if ~FacesTest(jj)
            stmp = s;
            stmp.vrtx(jj,:) = [];
            stmp = S1D(stmp);
            if stmp.nvrtx == 1
                d = norm(stmp.lambda(1)*stmp.vrtx(1,:));
            else
                d = norm(sum(stmp.lambda(1:stmp.nvrtx)'.*stmp.vrtx(1:stmp.nvrtx,:)));
            end
        if d < d_
            smin = stmp;
            d_ = d;
        end    
        end
    end
    s = smin;
end
end



function s = S1D(s)
% S1D
b = s.vrtx(1,:);
a = s.vrtx(2,:);
t = b - a;
nu_fabs = abs(t);

indexI = 2;
if nu_fabs(1) > nu_fabs(2)
    if nu_fabs(1) > nu_fabs(3)
        indexI = 1;
    else
        indexI = 3;
    end
elseif nu_fabs(1) < nu_fabs(2)
    if nu_fabs(2) > nu_fabs(3)
        indexI = 2;
    else
        indexI = 3;
    end
elseif nu_fabs(1) < nu_fabs(3)
    indexI = 3;
elseif nu_fabs(2) < nu_fabs(3)
    indexI = 3;
end

pt = (b*t')*(a(indexI) - b(indexI))/(t*t') + b(indexI);

det_ap = a(indexI) - pt;
det_pb = pt - b(indexI);

FacetsTest(1) = CompareSigns(t(indexI), - det_ap);
FacetsTest(2) = CompareSigns(t(indexI), - det_pb);

if sum(FacetsTest) == 2
    s.lambda(1) = - det_ap/t(indexI);
    s.lambda(2) = 1 - s.lambda(1);
%     s.wids = [1 2];
    s.nvrtx = 2;
elseif FacetsTest(1) == 0
    s.lambda = 1;
%     s.wids = 1;
    s.nvrtx = 1;
    s.vrtx(1,:) = s.vrtx(2,:);
else
    s.lambda = 1;
%     s.wids = 2;
    s.nvrtx = 1;
end
end 

function flag = CompareSigns(a,b)
    if a>0 && b>0
        flag = 1;
    elseif a<0 && b<0
        flag = 1;
    else
        flag = 0;
    end
end



function [x, distance] = findDistance(A)
    num = size(A,2);
    H = 2*(A'*A);
    Aeq = ones(1,num);
    beq = 1;
    lb = zeros(num,1);
    ub = ones(num,1);
    options = optimoptions('quadprog','Display','off');
    [x,fval] = quadprog(H,[],[],[],Aeq,beq,lb,ub,[],options);
    distance = sqrt(abs(fval));
end