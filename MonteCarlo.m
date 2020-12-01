function [msltn, rsdl2norm] = MonteCarlo(fwdfun, pararange, obsvstn, obsvdat)
%
% [msltn, rsdl2norm] = MonteCarlo(fwdfun, pararange, obsvstn, obsvdat)
%
% This is a Monte Carlo inversion program for the problem
% obsvdat = fwdfun(msltn, obsvstn) in the range restrain pararange of
% model parameters.
%
% Written by Tche L. at USTC, 2015/12.
%
% msltn: a vector whose size is [paranum, 1],
%        the model solution given by this program.
% rsdl2norm: a constant variable,
%            the 2-norm of the residual for the model msltn.
%
% fwdfun: a function handle,
%         the forward function.
% pararange: a matrix whose size is [paranum, 2],
%            every row is corresponding to a model parameter, the 1st column is
%            its minimum and the 2nd column is its maximum.
% obsvstn: a vector whose size is [datanum, ~],
%          the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum, 1],
%          the observe data used to inverse the problem
%          obsvdat = fwdfun(msltn, obsvstn).

[paranum, ~] = size(pararange);                                                % paranum is the number of model parameters.
itsmax = 100000;                                                               % the maximum time of Monte Carlo iteration.

modelset = NaN*ones(paranum, itsmax);                                          % the generated model set, every column represents a model.
r2normset = Inf*ones(1, itsmax);                                               % the set of the 2-norm of the residual for all models of the modelset.

for its = 1:1:itsmax
  randnum = rand(paranum, 1);                                                  % a random number vector in [0, 1].
  modelset(:, its) = pararange(:, 1) ...
    + randnum.*(pararange(:, 2) - pararange(:, 1));
  prdcdat = fwdfun(modelset(:, its), obsvstn);                                 % the predict data when using the itsth column of modelset as a model.
  r2normset(its) = norm(obsvdat - prdcdat, 2);
end

rsdl2norm = min(r2normset);
msltn = modelset(:, find(r2normset == rsdl2norm, 1, 'first'));

end
