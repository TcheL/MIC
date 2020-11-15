function [modelsltn, rsdl2norm] = HighorderGridsearch(fwdfun, pararange, ...
  parasmplnum, obsvstn, obsvdat)
%
% [modelsltn, rsdl2norm] = HighorderGridsearch(fwdfun, pararange, ...
%   parasmplnum, obsvstn, obsvdat)
%
% This is a high-order (multi-dimension) grid search inversion program for the
% problem obsvdat = fwdfun(modelsltn, obsvstn) in the range restrain pararange
% of model parameters.
% Note that this program can be used to a problem where the number of model
% parameters is arbitrary.
%
% Written by Tche L. at USTC, 2016/01.
%
% modelsltn: a vector whose size is [paranum, 1],
%            the model solution given by this program.
% rsdl2norm: a constant variable,
%            the 2-norm of the residual for the model modelsltn.
%
% fwdfun: a function handle,
%         the forward function.
% pararange: a matrix whose size is [paranum, 2],
%            every row is corresponding to a model parameter, the 1st column
%            is its minimum and the 2nd column is its maximum.
% parasmplnum: a vector whose size is [paranum, 1],
%              the number of the sampling points of every parameter.
% obsvstn: a vector whose size is [datanum, ~],
%          the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum, 1],
%          the observe data used to inverse the problem
%          obsvdat = fwdfun(modelsltn, obsvstn).

paranum = length(parasmplnum);                                                 % the number of the model parameters.

nodalpointnum = cumprod(parasmplnum)./parasmplnum;                             % the number of nodal points in the difference dimesion space; but note that the ith element represents total nodal points of the (i - 1)-dimesion space.
parasmplintv = (pararange(:, 2) - pararange(:, 1))./(parasmplnum - 1);         % the sample interval of every model parameter.

r2normset = Inf*ones(1, nodalpointnum(paranum));                               % the set of the 2-norm of the residual for all models.

for i = 1:prod(parasmplnum)
  smplpointpos = calpointpos(i, paranum, nodalpointnum);                       % the position coordinate of the ith calculate point in the multi-dimesion model space.
  model = pararange(:, 1) + parasmplintv.*(smplpointpos - 1);                  % the constructed model by using the position coordinate of the calculate point.
  prdcdat = fwdfun(model, obsvstn);                                            % the predict data when using the model.
  r2normset(i) = norm(obsvdat - prdcdat, 2);
end

bestpointordernum = find(r2normset == min(r2normset), 1, 'first');             % the order number (sequence number) of the best point where the residual is the least.
bestpointpos = calpointpos(bestpointordernum, paranum, nodalpointnum);         % the position coordinate of the best point in the multi-dimesion model space.

modelsltn = pararange(:, 1) + parasmplintv.*(bestpointpos - 1);
rsdl2norm = r2normset(bestpointordernum);

end

function pointposition = calpointpos(ordernum, paranum, nodalpointnum)
%
% pointposition = calpointpos(ordernum, paranum, nodalpointnum)
%
% This is a calculate point position coordinate program called by the
% multi-dimesion grid search program HighorderGridsearch.
%
% Written by Tche L. at USTC, 2016/01.

% pointposition: a vector whose size is [paranum, 1],
%                the point position coordinate calculated by this funtion.

% ordernum: a constant variable,
%           the order number (sequence number) of the calculate point.
% paranum: a constant variable,
%          the number of the model parameters,
%          i.e. paranum = length(nodalpointnum).
% nodalpointnum: a vector whose size is [paranum, 1],
%                the number of nodal points in the difference dimesion space;
%                but note that the ith element represents total nodal points of
%                the (i - 1)-dimesion space.

pointposition = NaN*ones(paranum, 1);

for i = paranum:-1:2
    if(mod(ordernum, nodalpointnum(i)) == 0)
        pointposition(i) = fix(ordernum/nodalpointnum(i));
        ordernum = nodalpointnum(i);
    else
        pointposition(i) = fix(ordernum/nodalpointnum(i)) + 1;
        ordernum = mod(ordernum, nodalpointnum(i));
    end
end

pointposition(1) = ordernum;

end
