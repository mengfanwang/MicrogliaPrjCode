% clc;clear;close all;
% 
% line_ind = [1 1 0; 2 1 0; 1 2 0; 2 2 0;];
% line_ind = line_ind';

function distance = gjk_v1(line_ind)
s.vrtx = zeros(4,3);
s.vrtx(1,:) = line_ind(:,1)';
vk = line_ind(:,1);
s.nvrtx = 1;
distance = norm(line_ind(:,1));

while vk'*vk >= 1e-4 && s.nvrtx < 4
wk_ind = findwk(line_ind,vk);
Wk = line_ind(:,wk_ind);

% [~, ~, distance] = subalgorithm(Y);
if ismember(Wk',s.vrtx(1:s.nvrtx,:),'row') ||  vk'*vk - vk'*Wk <= 1e-8 * (vk'*vk)
    break
end
s.vrtx(s.nvrtx+1,:) = line_ind(:,wk_ind)';
s.nvrtx = s.nvrtx + 1;

s = subalgorithm(s);


if s.nvrtx == 1
    vk = s.lambda.*s.vrtx(1:s.nvrtx,:)';
    distance = norm(s.lambda(1)*s.vrtx(1,:));
else
    vk = sum(s.lambda(1:s.nvrtx)'.*s.vrtx(1:s.nvrtx,:))';
    distance = norm(sum(s.lambda(1:s.nvrtx)'.*s.vrtx(1:s.nvrtx,:)));
end
end
end

% function distance = gjk_v1(line_ind)
% Y = 1;
% y = line_ind(:,1);
% vk = y;
% 
% while vk'*vk >= 1e-4
% wk_ind = findwk(line_ind,vk);
% [~, distance] = findDistance(line_ind(:,Y));
% if ismember(wk_ind,Y) || vk'*vk - vk'*line_ind(:,wk_ind) <= 1e-8 * (vk'*vk)
%     break
% end
% Y = [Y wk_ind];
% 
% % sub algorithm
% if length(Y) == 1
%     % 1-simplex
%     lambda = 1;
%     Wk = Y;
% else
%     % >2-simplex
%     [lambda, distance] = findDistance(line_ind(:,Y));
%     Wk = Y(lambda > 0);
% end
% vk = sum(lambda(lambda>0)'.*line_ind(:,Wk),2);
% Y = Wk;
% end
% end

function wk_ind = findwk(line_ind,vk)
temp_support = -inf;
for ii = 1:size(line_ind,2)
    if line_ind(:,ii)'*(-vk) > temp_support
        wk_ind = ii;
        temp_support = line_ind(:,ii)'*(-vk);
    end
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