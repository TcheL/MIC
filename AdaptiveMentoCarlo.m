function [msltn,rsdl2norm] = AdaptiveMentoCarlo(fwdfun,pararange,obsvstn,obsvdat)
% [msltn,rsdl2norm] = AdaptiveMentoCarlo(fwdfun,pararange,obsvstn,obsvdat)
% This is a adaptive Mento Carlo inversion program for the problem obsvdat = fwdfun(msltn,obsvstn) in the range restrain pararange of model parameters.
% Written by Tche.L. from USTC, 2015.12.
%
% msltn: a vector whose size is [paranum,1], the model solution given by this program.
% rsdl2norm: a constant variable, the 2-norm of the residual for the model msltn.
%
% fwdfun: a function handle, the forward function.
% pararange: a matrix whose size is [paranum,2], every row is corresponding to a model parameter, the 1st column is its minimum and the 2nd column is its maximum.
% obsvstn: a vector whose size is [datanum,~], the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum,1], the observe data used to inverse the problem obsvdat = fwdfun(msltn,obsvstn).

[paranum,~] = size(pararange);                                          % paranum is the number of model parameters.
modelnum = 1000;                                                        % the number of the generated model in every iteration step, i.e. do modelnum sets of parallel Mento Carlo search.
itsmax = 100;                                                           % the maximum time of adaptive Mento Carlo iteration.

modelbestprev = pararange(:,1)*ones(1,modelnum);                        % the previous best models, every column respresents the best model in all previous iteration steps for a set of parallel Mento Carlo search.
r2normsetprev = Inf*ones(1,modelnum);                                   % the set of the 2-norm of the residual for all models of the modelbestprev.
modelset = NaN*ones(paranum,modelnum);                                  % the new-generated model set, every column represents a model.
r2normset = Inf*ones(1,modelnum);                                       % the set of the 2-norm of the residual for all models of the modelset.

for its = 1:1:itsmax
    k = 2*its - 1;                                                      % the step-length attenuating exponent factor.
    for i = 1:1:modelnum
        while(1)
            randnum = rand(paranum,1);                                  % a random number vector in [0,1].
            modelset(:,i) = modelbestprev(:,i) + ...
                (2*randnum - 1).^k.*(pararange(:,2) - pararange(:,1));
            if(sum(modelset(:,i) >= pararange(:,1)) + ...
                    sum(modelset(:,i) <= pararange(:,2)) == 2*paranum)
                % if all model parameters are not less than their avaliable minimum and not greater than their avaliable maximum.
                break;
            end
        end
        prdcdat = fwdfun(modelset(:,i),obsvstn);                        % the predict data when using the ith column of modelset as a model.
        r2normset(i) = norm(obsvdat - prdcdat,2);
        if(r2normset(i) < r2normsetprev(i))
            modelbestprev(:,i) = modelset(:,i);
            r2normsetprev(i) = r2normset(i);
        end
    end
end

rsdl2norm = min(r2normsetprev);
msltn = modelbestprev(:,find(r2normsetprev == rsdl2norm,1,'first'));

end