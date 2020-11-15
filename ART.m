function [gridmodel, rsdl2norm] = ART(Gmatrix, traveltime, epsr2norm, epsupdate)
%
% [gridmodel, rsdl2norm] = ART(Gmatrix, traveltime, epsr2norm, epsupdate)
%
% This is a ART (Algebraic Reconstruction Technique) program that solve the
% tomographic reconstruction problem traveltime = Gmatrix*gridmodel.
%
% Written by Tche L. at USTC, 2015/12.
%
% gridmodel: a vector whose size is [paranum, 1],
%            the reconstructed grid model cell parameters.
% rsdl2norm: a constant variable,
%            the 2-norm of the residual of the model gridmodel.
%
% Gmatrix: a matrix whose size is [tnum, paranum],
%          the ray path matrix; Gmatrix(i, j) is l only when the ith ray
%          crosses the jth cell, otherwise it is 0; l is the length of the ith
%          ray in the jth cell.
% traveltime: a vector whose size is [tnum, 1],
%             the travel time data.
% epsr2norm: a small positive constant variable,
%            the residual's 2-norm to allow stop iterating.
% epsupdate: a small positive constant variable,
%            the maximum model update amount to allow stop iterating.
%
% Ref.: The Algorithm 6.2 on the P.147 of the textbook Parameter Estimation And
%       Inverse Problems (Second Edition. Richard C. Aster, Brian Borchers,
%       Clifford H. Thurber. 2011. Academic Press.)

[tnum, paranum] = size(Gmatrix);                                               % tnum is the number of the travel time data; paranum is the number of the model parameters.

itsmax = 10000;

updatecoeffi = (Gmatrix ~= 0);                                                 % the update coefficient matrix; updatecoeffi(i, j) is 1 when the jth cell in the ith ray path, otherwise it is 0.

its = 0;                                                                       % the iterating time.
raylen = sum(Gmatrix, 2);                                                      % the length of every ray path.
gridmodel = zeros(paranum, 1);

while(1)
  its = its + 1;
  maxmodelupdate = 0;                                                          % the maximum model update amounts.
  for i = 1:1:tnum
    for j = 1:1:paranum
      iprdctime = Gmatrix(i, :)*gridmodel;                                     % the ith predict travel time datum.
      modelupdate = (iprdctime - traveltime(i))*updatecoeffi(i, j) ...         % the model update amounts.
        /raylen(i);
      gridmodel(j) = gridmodel(j) - modelupdate;
      maxmodelupdate = max(maxmodelupdate, abs(modelupdate));
    end
  end
  rsdl2norm = norm(Gmatrix*gridmodel - traveltime, 2);
  if(rsdl2norm <= epsr2norm || maxmodelupdate <= epsupdate)
    break;
  elseif(its >= itsmax)
    warning('MATLAB:MICMaximumIterationTime', ['The time of ART ', ...
      'iterating comes to %d, and the ART iteration has been forced ', ...
      'to exit, as a result, the solution may be not enough accurate.'], ...
      its);
    break;
  end
end

end
