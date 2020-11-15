function [modelsltn, rsdl2norm] = LevenbergMarquardt(fwdfun, model, ...
  pararange, obsvstn, obsvdat, weightfactor, TheoreJacFun, ...
  steplenfinitediffer, damping, allowquitavemodification)
%
% [modelsltn, rsdl2norm] = LevenbergMarquardt(fwdfun, model, pararange, ...
%   obsvstn, obsvdat, weightfactor, TheoreJacFun, [], damping, ...
%   allowquitavemodification)
% [modelsltn, rsdl2norm] = LevenbergMarquardt(fwdfun, model, pararange, ...
%   obsvstn, obsvdat, weightfactor, [], steplenfinitediffer, damping, ...
%   allowquitavemodification)
%
% This is a Levenberg-Marquardt inversion program for the problem
% obsvdat = fwdfun(msltn, obsvstn).
%
% Written by Tche L. at USTC, 2016/01.
%
% modelsltn: a vector whose size is [paranum, 1],
%            the model solution given by this program.
% rsdl2norm: a constant variable,
%            the 2-norm of the residual for the model msltn.
%
% fwdfun: a function handle,
%         the forward function.
% model: a vector whose size is [paranum, 1],
%        the initial model.
% pararange: a matrix whose size is [paranum, 2],
%            every row is corresponding to a model parameter, the 1st column is
%            its minimum and the 2nd column is its maximum.
% obsvstn: a vector whose size is [datanum, ~],
%          the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum, 1],
%          the observe data used to inverse the problem
%          obsvdat = fwdfun(modelsltn,obsvstn).
% weightfactor: a vector whose size is [datanum, 1],
%               the weight factor of data.
% TheoreJacFun: a function handle,
%               the function that calculates theoretical Jacobian.
% steplenfinitediffer: a vector whose size is [paranum, 1],
%                      all of that are very small positive, the difference
%                      step-length of the finite-difference approximately
%                      calculating Jacobian.
% damping: a positive constant variable,
%          the initial damping factor.
% allowquitavemodification: a very small positive constant variable,
%                           the minimum amount of modification to allow to quit.
%
% Ref.: The equation (9.13) on the P.220 to the equation (9.30) on the P.222 of
%       the textbook Parameter Estimation And Inverse Problems (Second Edition.
%       Richard C. Aster, Brian Borchers, Clifford H. Thurber. 2011.
%       Academic Press.)

itsmax = 10000;                                                                % the maximum iterating time.
dampingmax = 1.0e+10;                                                          % the maximum value of the damping factor.
dampingmin = 1.0e-12;                                                          % the minimum value of the damping factor.
zoomfactor = 5;

paranum = length(model);                                                       % the number of model parameters.
prdcdatinit = fwdfun(model, obsvstn);                                          % the predict data of the initial model.
weightresinit = (prdcdatinit - obsvdat)./weightfactor;                         % the weight residual of prdcdatinit.

its = 0;                                                                       % the current iterating time.

while(1)
  its = its + 1;
  fprintf('its = %d, program is running ...\n', its);
  JacMat = CalculateJac(TheoreJacFun,steplenfinitediffer, ...                  % the Jacobian matrix.
    fwdfun, obsvstn, obsvdat, weightfactor, model);
  modelupdateamount = (JacMat'*JacMat + damping*eye(paranum)) ...
    \( - JacMat'*weightresinit);
  model = model + modelupdateamount;
  model(model < pararange(:, 1)) = pararange(model < pararange(:, 1), 1);      % If the model parameter is less than its lower bound, we let it equal to its lower bound value.
  model(model > pararange(:, 2)) = pararange(model > pararange(:, 2), 2);      % If the model parameter is lager than its upper bound, we let it equal to its upper bound value.
  prdcdatnew = fwdfun(model, obsvstn);                                         % the predict data of the new model.
  if(mean(abs(modelupdateamount)) <= allowquitavemodification)
    break;
  elseif(its >= itsmax)
    warning('MATLAB:MICMaximumIterationTime', ['The iterating time of ', ...
      'Levenberg-Marquardt iteration comes to %d, and force the loop to ', ...
      'exit. So the solution may be not enough accurate.'], its);
    break;
  end
  weightresnew = (prdcdatnew - obsvdat)./weightfactor;                         % the weight residual of prdcdatnew.
  convergeratio = norm(weightresnew, 2)/norm(weightresinit, 2);                % the ratio of the new weight residual to the initial one; it expresses the converging condition.
  if(convergeratio < 0.1)
    % if converge well
    if(damping < dampingmax)
      damping = damping*zoomfactor;
    else
      damping = dampingmax;
    end
  elseif(convergeratio > 1)
    % if not converge
    if(damping > dampingmin)
      damping = damping/zoomfactor;
    else
      damping = dampingmin;
    end
  end
  weightresinit = weightresnew;
end

modelsltn = model;
rsdl2norm = norm(prdcdatnew - obsvdat, 2);

end

function Jacobian = CalculateJac(TheoreJacFun, steplenfinitediffer, fwdfun, ...
  obsvstn, obsvdat, weightfactor, modelpoint)
%
% Jacobian = CalculateJac(TheoreJacFun, [], [], obsvstn, [], [], modelpoint)
% Jacobian = CalculateJac([], steplenfinitediffer, fwdfun, obsvstn, obsvdat, ...
%   weightfactor, modelpoint)
%
% This is a program used to calculate Jacobain by a Levenberg-Marquardt program.
%
% Written by Tche L. at USTC, 2016/01.

% Jacobian: a matrix whose size is [datanum, paranum],
%           the calculated Jacobian.

% TheoreJacFun: a function handle,
%               the function that calculates theoretical Jacobian.
% steplenfinitediffer: a vector whose size is [paranum, 1],i
%                      all of that are very small positive, the difference i
%                      step-length of the finite-difference approximately i
%                      calculating Jacobian.
% fwdfun: a function handle,
%         the forward function.
% obsvstn: a vector whose size is [datanum, ~],
%          the observe station point x for the data.
% obsvdat: a vector whose size is [datanum, 1],
%          the observe data used to inverse this problem.
% weightfactor: a vector whose size is [datanum, 1],
%               the weight factor of data.
% modelpoint: a vector whose size is [paranum, 1],
%             the model point where we calculate the Jacobian;
%             Jacobain(modelpoint), i.e. J(m).

if( ~ isempty(TheoreJacFun))
  Jacobian = TheoreJacFun(modelpoint, obsvstn);
  % The above program statement can be properly modified based on
  % the actual situation by user.
elseif( ~ isempty(steplenfinitediffer))
  paranum = length(modelpoint);                                                % the number of model parameters.
  datacalpoint = fwdfun(modelpoint, obsvstn);                                  % the forward data of the calculated point.
  weirescalpoint = (datacalpoint - obsvdat)./weightfactor;                     % the weight residual of datacalpoint.
  Jacobian = NaN*ones(length(obsvdat), paranum);
  for i = 1:1:paranum
    modelnearbypoint = modelpoint;                                             % the model point nearby the calculated point.
    modelnearbypoint(i) = modelpoint(i) + steplenfinitediffer(i);
    datanearbypoint = fwdfun(modelnearbypoint, obsvstn);                       % the forward data of modelnearbypoint
    weiresnearbypoint = (datanearbypoint - obsvdat)./weightfactor;             % the weight residual of datanearbypoint.
    Jacobian(:,i) = (weiresnearbypoint - weirescalpoint) ...
      ./steplenfinitediffer(i);
  end
else
  error(['The function CalculateJac goes wrong because both the two ', ...
    'parameters TheoreJacFun and steplenfinitediffer do not exist!');
end

end
