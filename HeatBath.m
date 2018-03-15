function [msltn,rsdl2norm] = HeatBath(fwdfun,Tempseries,pararange,parasmplnum,obsvstn,obsvdat)
% [msltn,rsdl2norm] = HeatBath(fwdfun,Tempseries,pararange,parasmplnum,obsvstn,obsvdat)
% This is a Heat Bath inversion program for the problem obsvdat = fwdfun(msltn,obsvstn) in the range restrain pararange of model parameters.
% Written by Tche.L. from USTC, 2015.12.
%
% msltn: a vector whose size is [paranum,1], the model solution given by this program.
% rsdl2norm: a constant variable, the 2-norm of the residual for the model msltn.
%
% fwdfun: a function handle, the forward function.
% Tempseries: a vector whose size is [Tempnum,1], the annealing temperature series.
% pararange: a matrix whose size is [paranum,2], every row is corresponding to a model parameter, the 1st column is its minimum and the 2nd column is its maximum.
% parasmplnum: a vector whose size is [paranum,1], the sample number of every model parameter in its range.
% obsvstn: a vector whose size is [datanum,~], the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum,1], the observe data used to inverse the problem obsvdat = fwdfun(msltn,obsvstn).
%
% Ref.: Heat Bath Algorithm on the P.92 of the textbook Global Optimization Methods in Geophysical Inversion (Second Edition. Mrinal K. Sen, Paul L. Stoffa. 2012. Cambridge University Press.)

paranum = length(parasmplnum);                                                      % the number of model parameters.
Tempnum = length(Tempseries);                                                       % the number of temperature sample points.
randmovenum = randi([50 100],1,Tempnum);                                            % the number of random moves of every temperature sample point.
randmovesum = sum(randmovenum);                                                     % the sum of randmovenum.
parasmplintv = (pararange(:,2) - pararange(:,1))./(parasmplnum - 1);                % the sample interval of every model parameter.
parasmplmax = max(parasmplnum);                                                     % the maximum value of parasmplnum.

modelinit = pararange(:,1) + ...                                                    % the initial model.
    rand(paranum,1).*(pararange(:,2) - pararange(:,1));
modelset = NaN*ones(paranum,randmovesum);                                           % the model set generated in every iteration, every column represents a model.
r2normset = Inf*ones(1,randmovesum);                                                % the set of the 2-norm of the residual for all models of the modelset.

modelparapoint = NaN*ones(paranum,parasmplmax);                                     % the available sample points of every model parameter.
for k = 1:1:paranum
    modelparapoint(k,:) = pararange(k,1) + (0:1:parasmplnum(k) - 1).*parasmplintv(k);
end
errorenergy = Inf*ones(1,parasmplmax);                                              % the error energy of a certain model parameter, i.e. the 2-norm numbers of the residual when using a series of available parameter sample points for a certain model parameter.
Gibbsprob = zeros(1,parasmplmax);                                                   % the Gibbs probabilities of a series of available sample points for a certain model parameter.
startpointer = cumsum(randmovenum) - randmovenum;                                   % the position pointer of the starting point when do temperature loop.

for i = 1:1:Tempnum
    for j = 1:1:randmovenum(i)
        for k = 1:1:paranum
            for l = 1:1:parasmplnum(k)
                modelinit(k) = modelparapoint(k,l);
                prdcdat = fwdfun(modelinit,obsvstn);                                % the predict data when using modelinit as a model.
                errorenergy(l) = norm(obsvdat - prdcdat,2);
                Gibbsprob(l) = exp( - errorenergy(l)/Tempseries(i));
            end
            Gibbsprobsum = sum(Gibbsprob(1:parasmplnum(k)));                        % the sum of the Gibbs probabilities for a certain model parameter.
            if(Gibbsprobsum == 0)
                % if each errorenergy is very large, the sum of Gibbsprob is close to 0, the latter procedure will go wrong. Thus, if so, immediately jump out to do next loop.
                continue;
            end
            Gibbsprobcum = cumsum(Gibbsprob(1:parasmplnum(k))./Gibbsprobsum);       % the cumulative Gibbs probability for a certain model parameter.
            iprob = find(rand <= Gibbsprobcum,1,'first');                           % the iprobth probability interval is selected, i.e. the iprobth model parameter sample point of modelparapoint is selected for a certain model parameter.
            modelinit(k) = modelparapoint(k,iprob);
        end
        modelset(:,startpointer(i) + j) = modelinit;
        r2normset(startpointer(i) + j) = errorenergy(iprob);
    end
end

rsdl2norm = min(r2normset);
msltn = modelset(:,find(r2normset == rsdl2norm,1,'first'));

end