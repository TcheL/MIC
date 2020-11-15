function [model, rsdl2norm] = Kaczmarz(Gmatrix, data, epsr2norm, epsupdate)
%
% [model, rsdl2norm] = Kaczmarz(Gmatrix, data, epsr2norm, epsupdate)
%
% This is a Kaczmarz's Algorithm program that solve the problem
% data = Gmatrix*model.
%
% Written by Tche L. at USTC, 2015/12.
%
% model: a vector whose size is [paranum, 1],
%        the model parameter solutions given by this program.
% rsdl2norm: a constant variable,
%            the 2-norm of the residual of the model.
%
% Gmatrix: a matrix whose size is [datanum, paranum],
%          the matrix G for d = G*m.
% data: a vector whose size is [datanum, 1],
%       the data for d = G*m.
% epsr2norm: a small positive constant variable,
%            the residual's 2-norm to allow stop iterating.
% epsupdate: a small positive constant variable,
%            the average model update amount to allow stop iterating.
%
% Ref.: The Algorithm 6.1 on the P.146 of the textbook Parameter Estimation
%       And Inverse Problems (Second Edition. Richard C. Aster, Brian Borchers,
%       Clifford H. Thurber. 2011. Academic Press.)

[datanum, paranum] = size(Gmatrix);                                            % datanum is the number of data, paranum is the number of the model parameters.

itsmax = 10000;                                                                % the maximum time of iterating.

its = 0;                                                                       % the current time of iterating.
model = zeros(paranum, 1);
modelupdate = zeros(paranum, datanum);                                         % the model update amounts.

while(1)
  its = its + 1;
  for i = 1:1:datanum
    modelupdate(:, i) = (Gmatrix(i, :)*model - data(i))/(Gmatrix(i, :) ...
      *Gmatrix(i, :)')*Gmatrix(i, :)';
    model = model - modelupdate(:, i);
  end
  rsdl2norm = norm(Gmatrix*model - data, 2);
  if(rsdl2norm <= epsr2norm || sum(abs(modelupdate(:)))/paranum <= epsupdate)
    break;
  elseif(its >= itsmax)
    warning('MATLAB:MICMaximumIterationTime', ['The time of Kaczmarz''s ', ...
      'iterating comes to %d, and the Kaczmarz''s iteration has been ', ...
      'forced to exit, as a result, the solution may be not enough ', ...
      'accurate.'], its);
    break;
  end
end

end
