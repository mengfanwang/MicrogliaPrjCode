% clc;clear;close all;dbstop if error
function [stabilizeFunction, variance] = convexOptimization(his,parameters)

    % This function is used to stabilize the variance from the given
    % histograms.
    %
    % function stabilize_function = ConvexOptimization(his,parameters)
    %
    % Input: his ------------ An estimated n-by-m histogram matrix. n is 
    %                         the number of histograms and m is the number 
    %                         of bins in each histogram. Each row is a 
    %                         histogram whose real signals are the
    %                         corresponding entry in histCenters.
    %
    %        parameters: A necessary structure contains some useful 
    %                    information.
    %        .var_ratio ----- A n dimension vector to control the weight of
    %                         the variance. This function will make 
    %                         variance./var_ratio be almost the same. The
    %                         default value is all ones.
    %        .histCenters --- A n dimension vector shows the center or the
    %                         real signal of every histogram.
    %        .binCenters ---- A m dimension vector shows the center or the
    %                         signal of every bin in the histograms.
    %        .step ---------- A number>0 to decide use the optimzation
    %                         steps. If .step = 0, the algorithm will only
    %                         use the first step optimization. Otherwise it
    %                         will use the second step and circulate .step
    %                         times.
    %                        
    %        .algorithm 
    % 
    % Output: stabilize_function
    %               --------- A m dimension vector used to transform the 
    %                         signal. Every entry is the transformed value
    %                         of binCenters. The first and the last one is
    %                         the same as those in binCenters.
    



% zero_ind = find(sum(his)==0);
% parameters.binCenters(zero_ind) = [];
% his(:,zero_ind) = [];
nargin = 2;

if any(sum(his,2) < 1 - 1e-6) || any(sum(his,2) > 1 + 1e-6)
    warning('The sum of some histograms may be not 1.');
    his = his./sum(his,2);
end
if nargin ~= 2
 error('The number of inputs is incorrect.');
end
if ~isfield(parameters,'histCenters')
error('parameters.histCenters is necessary.');
end
if ~isfield(parameters,'binCenters')
error('parameters.binCenters is necessary.');
end
histCenters = parameters.histCenters;
binCenters = parameters.binCenters; 
if min(histCenters) < min(binCenters) || max(histCenters) > max(binCenters)
    error('The range of histCenters should not exceed the range of binCenters.');
end
numOfHists = length(histCenters);
numOfBins = length(binCenters);
if ~isfield(parameters,'var_ratio')
    var_ratio = ones(1,numOfHists);
else
    var_ratio = parameters.var_ratio;
end
if ~isfield(parameters,'step')
    step = 1;
else
    step = parameters.step;
end
if ~isfield(parameters,'algorithm')
    algorithm = 'aaaaa';
else
    algorithm = 'bbbbb';
end

% optimization part first step
tic;
disp('Start First Step...');

% MOSEK Part`
addpath('C:\Program Files\Mosek\9.2\toolbox\R2015a');
Hess = cell(1,numOfHists);
for ii = 1:numOfHists
    Hess{ii} = zeros(numOfBins-1,numOfBins-1);
    left_ind = find(binCenters<histCenters(ii),1,'last'); % also the split variable
    right_ind = left_ind+1; 
    alpha = (histCenters(ii) - binCenters(left_ind))/(binCenters(right_ind) - binCenters(left_ind));
    beta = (binCenters(right_ind) - histCenters(ii))/(binCenters(right_ind) - binCenters(left_ind));
    if isempty(left_ind)
        left_ind = 0; right_ind = 1; alpha = 1; beta = 0;
    end
    for jj = 1:numOfBins
        if jj < left_ind
            Hess{ii}(jj:left_ind-1,jj:left_ind-1) = Hess{ii}(jj:left_ind-1,jj:left_ind-1) + his(ii,jj);
            Hess{ii}(left_ind,jj:left_ind-1) = Hess{ii}(left_ind,jj:left_ind-1) + alpha*his(ii,jj);
            Hess{ii}(jj:left_ind-1,left_ind) = Hess{ii}(jj:left_ind-1,left_ind) + alpha*his(ii,jj);
            Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + alpha^2*his(ii,jj);
        elseif jj == left_ind
            Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + alpha^2*his(ii,jj);
        elseif jj == right_ind
            if left_ind > 0
                Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + beta^2*his(ii,jj);
            end
        elseif jj > right_ind
            Hess{ii}(right_ind:jj-1,right_ind:jj-1) = Hess{ii}(right_ind:jj-1,right_ind:jj-1) + his(ii,jj);
            if left_ind > 0 
                Hess{ii}(left_ind,right_ind:jj-1) = Hess{ii}(left_ind,right_ind:jj-1) + beta*his(ii,jj);
                Hess{ii}(right_ind:jj-1,left_ind) = Hess{ii}(right_ind:jj-1,left_ind) + beta*his(ii,jj);
                Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + beta^2*his(ii,jj);
            end
        end
    end
    Hess{ii} = tril(Hess{ii})*2/var_ratio(ii);
    Hess{ii}(Hess{ii}<1e-16) = 0;
end
prob.qcsubk = [];
prob.qcsubi = [];
prob.qcsubj = [];
prob.qcval = [];
for ii = 1:numOfHists
    [row,column] = find(Hess{ii});
    prob.qcsubk = [prob.qcsubk; ii*ones(length(row),1)];
    prob.qcsubi = [prob.qcsubi; row];
    prob.qcsubj = [prob.qcsubj; column];
    prob.qcval = [prob.qcval; Hess{ii}(sub2ind(size(Hess{ii}),row,column))];
end
prob.c = [zeros(1,numOfBins-1) 1];
prob.buc = [zeros(numOfHists,1); binCenters(end)-binCenters(1)];
prob.blc = [-inf*ones(numOfHists,1); binCenters(end)-binCenters(1)];
prob.blx = zeros(numOfBins,1);
prob.a = repmat([zeros(1,numOfBins-1) -1],[numOfHists 1]);
prob.a = [prob.a; [ones(1,numOfBins-1) 0]];
param.MSK_IPAR_INTPNT_SOLVE_FORM = 'MSK_SOLVE_DUAL';
[~,res] = mosekopt('minimize',prob,param);

x = res.sol.itr.xx(1:end-1)';
stabilizeFunction = cumsum([binCenters(1) x]);
f = @(x)objectiveFunction(x,his,var_ratio,numOfHists,binCenters,histCenters);
variance = f(x);
plot(variance);
end

function F = objectiveFunction(x,his,var_ratio,numOfHists,binCenters,histCenters)
    for ii = 1:numOfHists
        left_ind = find(binCenters<histCenters(ii),1,'last');
        if isempty(left_ind)
            F(ii) = sum(his(ii,:).*([0 cumsum(x)]).^2)/var_ratio(ii);
        else
            beta = (binCenters(left_ind+1) - histCenters(ii))/(binCenters(left_ind+1) - binCenters(left_ind));
            F(ii) = sum(his(ii,:).*([0 cumsum(x)]-(sum(x(1:left_ind))-beta*x(left_ind))).^2)/var_ratio(ii);
        end
    end
end
