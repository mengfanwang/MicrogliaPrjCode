function [mu, sigma] = conventionalIntegralv1(interSectV, outerBorderV)

%%%%%%%%%%%%
% interSectV: the values of the pixels in the gap
% outerBorderV: the values of the pixels in the neighbor of the gap 
% groupid is the vector containing the group information for the ordered
% variables, with 1 for group 1 and -1 for group 2.
% M is the number of elements in group1 and N is the number of
% elements in the group 2

M = length(outerBorderV);
N = length(interSectV);
valueAndlabel = zeros(M+N,2);
n = M+N;
valueAndlabel(1:N,1) = interSectV;
valueAndlabel(N+1:n,1) = outerBorderV;
valueAndlabel(1:N,2) = -1;
valueAndlabel(N+1:n,2) = 1;

valueAndlabel = sortrows(valueAndlabel, 1);
weights = zeros(1,M+N);
weights(valueAndlabel(:,2) <0) = -(M+N)/N;
weights(valueAndlabel(:,2) >0) = (M+N)/M;

u = [1/(n+1):1/(n+1):n/(n+1)];
v = u;
J = weights;
Fmu = J.*norminv(u);
mu = 1/2*(1/(n+1))*(2*sum(Fmu) - Fmu(1) - Fmu(end));


F1 = J.*u./normpdf(norminv(u)); % Here the normal distribution is used
F2 = J.*(1-v)./normpdf(norminv(v));
F_matrix = F1' * F2;
F_matrix_1 = triu(F_matrix,1);
S1 = 1/6*(1/(n+1))^2*(2*sum(diag(F_matrix_1(1:end-1,2:end)))+ sum(diag(F_matrix_1(1:end-2,3:end)))...
    -F_matrix_1(1,2) - F_matrix_1(end-1,end)) + 1/4*(1/(n+1)^2)*(sum(diag(F_matrix_1(1:end-1,2:end)))...
    -F_matrix_1(1,2) - F_matrix_1(end-1,end) + 3*sum(diag(F_matrix_1(1:end-2,3:end))) - 2*F_matrix_1(1,3)...
    -2*F_matrix_1(end-2,end) + 4*sum(sum(F_matrix_1(2:end-4,5:end-1))) + 2*(sum(F_matrix_1(1,4:end-1))+sum(F_matrix_1(2:end-3,end)))...
    +F_matrix_1(1,end));

F_matrix_2 = triu(F_matrix,0);
S2 = 1/6*(1/(n+1))^2*(2*sum(diag(F_matrix_2)) + sum(diag(F_matrix_2(1:end-1,2:end)))...
    -F_matrix_2(1,1) - F_matrix_2(end,end)) + 1/4*(1/(n+1)^2)*(sum(diag(F_matrix_2)) - F_matrix_2(1,1) - F_matrix_2(end,end)...
    + 3*(sum(diag(F_matrix_2(1:end-1,2:end)))) - 2*F_matrix_2(1,2) - 2*F_matrix_2(end-1,end)...
    +4*sum(sum(F_matrix_2(2:end-3,4:end-1))) + 2*(sum(F_matrix_2(1,3:end-1))+ sum(F_matrix_2(2:end-2,end))));

% F_matrix_3 = diag(F_matrix);
% S1 = 1/4*(1/(n+1))^2*(F_matrix_3(1,1) + F_matrix_3(end,end) + F_matrix_3(1,end) + F_matrix_3(end,1)...
%     + 2*sum(F_matrix_3(2:end-1,1)) +2*sum(F_matrix_3(2:end-1,end)) + 2*sum(F_matrix_3(1,2:end-1)) ...
%     + 2*sum(F_matrix_3(end,2:end-1)) + 4*sum(sum(F_matrix_3(2:end-1,2:end-1))));

sigma = sqrt(max(S1+S2,0))/sqrt(n);
end

