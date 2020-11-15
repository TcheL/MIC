function example1()

%%
x = linspace(1, 10, 30)';
mtrue = [2.5, 1.5, 3]';
d = myfun(mtrue, x);

%% Levenberg-Marquardt Method
% minit = [1, 1, 1]';
% mrange = [ 0, 5;
%           -1, 4;
%            0, 6];
% weight = ones(length(x), 1);
% stepdiff = 0.0001.*ones(length(mtrue), 1);
% 
% [mlmq, rlmq] = LevenbergMarquardt(@myfun, minit, mrange, x, d, weight, ...
%   [], stepdiff, 1.0, 1.0e-6)

%% Grid Search Method
% mrange = [ 0, 5;
%           -1, 4;
%            0, 6];
% nsample = [31, 31, 31]';
% 
% [mhgs, rhgs] = HighorderGridsearch(@myfun, mrange, nsample, x, d)

%% Mento Carlo Method
% mrange = [ 0, 5;
%           -1, 4;
%            0, 6];
% nsample = [31, 31, 31]';
% 
% [mmtc, rmtc] = MentoCarlo(@myfun, mrange, x, d)
% [mamc, ramc] = AdaptiveMentoCarlo(@myfun, mrange, x, d)

%% Simulate Annealing Method
% temp = 100.0*exp(-0.1*(1:1:100));
% mrange = [ 0, 5;
%           -1, 4;
%            0, 6];
% nsample = [31, 31, 31]';
% minit = [0, 0, 0]';
% 
% [mhtb, rhtb] = HeatBath(@myfun, temp, mrange, nsample, x, d)
% [mmps, rmps] = Metropolis(@myfun, minit, temp, mrange, x, d)

%% Genetic Algorithm
% mrange = [ 0, 5;
%           -1, 4;
%            0, 6];
% genelen = [10, 10, 10]';
% maxfit = 100.0;

% [msga, rsga] = StandardGeneticAlgorithm(@myfun, mrange, genelen, x, d, ...
%   maxfit, 0.99, 0.95)
% [mmga, rmga] = ModifiedGeneticAlgorithm(@myfun, mrange, genelen, x, d, ...
%   maxfit, 0.99, 0.95)

end

%% forward function
function y = myfun(m, x)
% y = a*x^2 + b*x + c.

a = m(1);
b = m(2);
c = m(3);

y = x.^2.*a + x.*b + c;

end
