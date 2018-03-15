function [modelsltn,rsdl2norm] = ConjugateGradient(Gmatrix,data,epsr2norm,epsupdate)
% [model,rsdl2norm] = ConjugateGradient(Gmatrix,data,epsr2norm,epsupdate)
% This is a CG (Conjugate Gradient Method) program that solve the problem d = Gm, i.e. data = Gmatrix*model, that don't have to be a symmetric positive definite system of equations.
% Note that this program can similarly use CGLS (conjugate gradient least squares method) to solve the least squares problem min ||Gm - d||_2 for d = Gm.
% Written by Tche.L. from USTC, 2015,12.
%
% modelsltn: a vector whose size is [paranum,1], the model solutions given by this program.
% rsdl2norm: a constant variable, the 2-norm of the residual of the model modelsltn.
%
% Gmatrix: a matrix whose size is [datanum,paranum], the matrix G for d = G*m.
% data: a vector whose size is [datanum,1], the data for d = G*m.
% epsr2norm: a small positive constant variable, the residual's 2-norm to allow stop iterating.
% epsupdate: a small positive constant variable, the average model update amount to allow stop iterating.
%
% Ref.: The Algorithm 6.4 on the P.154 and the Algorithm 6.5 on the P.155 of the textbook Parameter Estimation And Inverse Problems (Second Edition. Richard C. Aster, Brian Borchers, Clifford H. Thurber. 2011. Academic Press.)

itsmax = 10000;                                                             % the maximum time of iterating.

if(~ IsSymPosDef(Gmatrix))
    % if the matrix Gmatrix is not a symmetric positive definite matrix.
    data = Gmatrix'*data;
    Gmatrix = Gmatrix'*Gmatrix;
end

[datanum,paranum] = size(Gmatrix);                                          % the number of data

its = 0;                                                                    % the current time of iterating.
betastepfactor = 0;                                                         % beta: the step length facotr.
pkbasisvector = zeros(datanum,1);                                           % pk: the basis vector of the x model.
xkmodel = zeros(paranum,1);                                                 % xk: the solution vector.
rklastresidual = data - Gmatrix*xkmodel;                                    % rk: the last residual vector.

while(1)
    its = its + 1;
    pkbasisvector = rklastresidual + betastepfactor*pkbasisvector;
    alphacoeffi = (rklastresidual'*rklastresidual) ...                      % alpha: the coefficient of these basis vector consist of the solution x.
        /(pkbasisvector'*Gmatrix*pkbasisvector);
    modelupdate = alphacoeffi*pkbasisvector;
    xkmodel = xkmodel + modelupdate;
    rknewresidual = rklastresidual - alphacoeffi*Gmatrix*pkbasisvector;
    if(norm(rknewresidual,2) <= epsr2norm || mean(abs(modelupdate)) <= epsupdate)
        break;
    elseif(its >= itsmax)
        warning('MATLAB:MyMaximumIterationTime','The time of Conjugate Gradient iterating comes to %d, and the Conjugate Gradient  iteration has been forced to exit, as a result, the solution may be not enough accurate.',its);
        break;
    end
    betastepfactor = (rknewresidual'*rknewresidual)/(rklastresidual'*rklastresidual);
    rklastresidual = rknewresidual;
end

modelsltn = xkmodel;
rsdl2norm = norm(Gmatrix*modelsltn - data,2);

end

function SpdIndex = IsSymPosDef(A)
% SpdIndex = IsSymPosDef(A)
% This is a program of judging a matrix whether it is a symmetric positive definite matrix.
% If A is a symmetric positive definite matrix, the function will return SpdIndex = 1, otherwise, SpdIndex = 0.
% Written by Tche.L. from USTC, 2015,12.

% SpdIndex: a constant variable, the symmetric-positive-definite index.

% A: a matrix, the matrix we want to judge whether it is a symmetric positive definite matrix.

[nrow,ncol] = size(A);                                                      % nrow is the row number of A, ncol is the column number of A.

if(nrow == ncol)
    eigenvalue = eig(A);                                                    % the eigenvalue of the matrix A.
    if(sum(sum(A' == A)) + sum(eigenvalue > 0) == nrow*ncol + nrow)
        SpdIndex = 1;
    else
        SpdIndex = 0;
    end
else
    SpdIndex = 0;
end

end
