function [gridmodel,rsdl2norm] = SIRT(Gmatrix,traveltime,epsr2norm,epsupdate)
% [gridmodel,rsdl2norm] = SIRT(Gmatrix,traveltime,epsr2norm,epsupdate)
% This is a SIRT (Simultaneous Iterative Reconstruction Technique) program that solve the tomographic reconstruction problem traveltime = Gmatrix*gridmodel.
% Written by Tche.L. from USTC, 2015,12.
%
% gridmodel: a vector whose size is [paranum,1], the reconstructed grid model cell parameters.
% rsdl2norm: a constant variable, the 2-norm of the residual of the model gridmodel.
%
% Gmatrix: a matrix whose size is [tnum,paranum], the ray path matrix; Gmatrix(i,j) is l only when the ith ray crosses the jth cell, otherwise it is 0; l is the length of the ith ray in the jth cell.
% traveltime: a vector whose size is [tnum,1], the travel time data.
% epsr2norm: a small positive constant variable, the residual's 2-norm to allow stop iterating.
% epsupdate: a small positive constant variable, the average model update amount to allow stop iterating.
%
% Ref.: The Algorithm 6.3 on the P.148 of the textbook Parameter Estimation And Inverse Problems (Second Edition. Richard C. Aster, Brian Borchers, Clifford H. Thurber. 2011. Academic Press.)

[~,paranum] = size(Gmatrix);                                            % tnum is the number of the travel time data; paranum is the number of the model parameters.

itsmax = 10000;                                                         % the maximum time of iterating.

updatecoeffi = (Gmatrix ~= 0);                                          % the update coefficient matrix; updatecoeffi(i,j) is 1 when the jth cell in the ith ray path, otherwise it is 0.

its = 0;                                                                % the iterating time.
raynum = sum(Gmatrix ~= 0,1)';                                           % the number of ray paths that pass through every model cell.
raylen = sum(Gmatrix,2);                                                % the length of every ray path.
cellnum = sum(Gmatrix ~= 0,2);                                          % the number of cells traversed by every ray path.
gridmodel = zeros(paranum,1);

while(1)
    its = its + 1;
    modelupdate = zeros(paranum,1);                                     % the model update amounts.
    for i = 1:1:size(Gmatrix,1)
        iprdctime = (Gmatrix(i,:) ~= 0)*gridmodel;                      % the ith predict travel time divided by the cell dimension.
        modelupdate = modelupdate + (traveltime(i)/raylen(i) - iprdctime/cellnum(i))*updatecoeffi(i, :)';
    end
    gridmodel = gridmodel + modelupdate./raynum;
    rsdl2norm = norm(Gmatrix*gridmodel - traveltime,2);
    if(rsdl2norm <= epsr2norm || mean(abs(modelupdate)) <= epsupdate)
        break;
    elseif(its >= itsmax)
        warning('MATLAB:MyMaximumIterationTime','The time of SIRT iterating comes to %d, and the SIRT iteration has been forced to exit, as a result, the solution may be not enough accurate.',its);
        break;
    end
end


end
