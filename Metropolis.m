function [msltn, rsdl2norm] = Metropolis(fwdfun, modelinit, Tempseries, ...
  pararange, obsvstn, obsvdat)
%
% [msltn, rsdl2norm] = Metropolis(fwdfun, modelinit, Tempseries, pararange, ...
%   obsvstn, obsvdat)
%
% This is a Metropolis inversion program for the problem
% obsvdat = fwdfun(msltn, obsvstn) in the range restrain pararange of
% model parameters.
%
% Written by Tche L. at USTC, 2015/12.
%
% msltn: a vector whose size is [paranum, 1],
%        the model solution given by this program.
% rsdl2norm: a constant variable,
%             the 2-norm of the residual for the model msltn.
%
% fwdfun: a function handle,
%         the forward function.
% modelinit: a vector whose size is [paranum, 1],
%            the initial model.
% Tempseries: a vector whose size is [Tempnum, 1],
%             the annealing temperature series.
% pararange: a matrix whose size is [paranum, 2],
%            every row is corresponding to a model parameter, the 1st column is
%            its minimum and the 2nd column is its maximum.
% obsvstn: a vector whose size is [datanum, ~],
%          the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum, 1],
%          the observe data used to inverse the problem
%          obsvdat = fwdfun(msltn, obsvstn).
%
% Ref.: Metropolis Algorithm on the P.86 of the textbook Global Optimization
%       Methods in Geophysical Inversion (Second Edition. Mrinal K. Sen,
%       Paul L. Stoffa. 2012. Cambridge University Press.)

paranum = length(modelinit);                                                   % the number of model parameters.
Tempnum = length(Tempseries);                                                  % the number of temperature sample points.
randmovenum = randi([50 100], 1, Tempnum);                                     % the number of random moves of every temperature sample point.

prdcdat = fwdfun(modelinit, obsvstn);                                          % the predict data when using the initial model as a model.
r2norminit = norm(obsvdat - prdcdat, 2);                                       % the 2-norm of the residual for the initial model.

for i = 1:1:Tempnum
  for j = 1:1:randmovenum(i)
    while(1)
      stepcoeffi = rand(paranum, 1) - 0.5;                                     % the random step coefficient, a random in [-0.5, 0.5].
      steplength = (pararange(:, 2) - pararange(:, 1)) ...                     % the random step length, less than half of the range span.
        .*rand(paranum, 1)./2;
      modelnew = modelinit + stepcoeffi.*steplength;                           % the new model.
      if(sum(modelnew >= pararange(:, 1)) ...
        + sum(modelnew <= pararange(:, 2)) == 2*paranum)
        % if all model parameters of modelnew are not less than their avaliable
        % minimum and not greater than their avaliable maximum
        break;
      end
    end
    prdcdat = fwdfun(modelnew, obsvstn);                                       % the predict data when using the new model as a model.
    r2normnew = norm(obsvdat - prdcdat, 2);                                    % the 2-norm of the residual for the new model.
    r2normdiff = r2normnew - r2norminit;                                       % the difference of the 2-norms of the two residuals for modelinit and modelnew.
    if(r2normdiff <= 0)
      modelinit = modelnew;
      r2norminit = r2normnew;
    else
      Gibbsprob = exp( - r2normdiff/Tempseries(i));                            % the Gibbs probability of the new model.
      if(rand < Gibbsprob)
        modelinit = modelnew;
        r2norminit = r2normnew;
      end
    end
  end
end

msltn = modelinit;
rsdl2norm = r2norminit;

end
