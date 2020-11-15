function [msltn, rsdl2norm] = ModifiedGeneticAlgorithm(fwdfun, pararange, ...
  paragenelen, obsvstn, obsvdat, fitnessmaxlimit, allowquitfitmax, ...
  allowquitfitavemaxratio)
%
% [msltn, rsdl2norm] = ModifiedGeneticAlgorithm(fwdfun, pararange, ...
%   paragenelen, obsvstn, obsvdat, fitnessmaxlimit, allowquitfitmax, ...
%   allowquitfitavemaxratio)
%
% This is a modified genetic algorithm inversion program for the problem
% obsvdat = fwdfun(msltn, obsvstn) in the range restrain pararange of
% model parameters.
% This program combines elements of SA into a new GA.
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
% paragenelen: a vector whose size is [paranum, 1],
%              the length of the gene fragment of every model parameter.
% obsvstn: a vector whose size is [datanum, ~],
%          the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum, 1],
%          the observe data used to inverse the problem
%          obsvdat = fwdfun(msltn, obsvstn).
% fitnessmaxlimit: a constant variable,
%                  the maximum value of fitness that can be tolerated.
% allowquitfitmax: a constant variable who is in [0, 1],
%                  the fitness value when the maximum fitness of a generation
%                  population is greater than it we can quit the genetic
%                  heredity iteration.
% allowquitfitavemaxratio: a constant variable who is in [0, 1],
%                          the ratio value when the ratio of the average fitness
%                          to the maximum fitness of a generation population is
%                          greater than it we can quit the genetic heredity
%                          iteration.
%
% Ref.: Genetic Algorithm on the P.135 of the textbook Global Optimization
%       Methods in Geophysical Inversion (Second Edition. Mrinal K. Sen,
%       Paul L. Stoffa. 2012. Cambridge University Press.)

paranum = length(paragenelen);                                                 % the number of model parameters.
popunum = 2^10;                                                                % the population number of every generation.
probcrossover = 0.7;                                                           % the probability of crossover.
probmutation = 0.05;                                                           % the probability of mutation.
ignmax = 100;                                                                  % the generation number of genetic heredity.
Tempupdate = 100;                                                              % the updating simulated annealing temperature.
Tempselect = 100;                                                              % the selecting simulated annealing temperature.

% Encoding (Coding of model parameters)
parasmplintv = (pararange(:, 2) - pararange(:, 1)) ...                         % the sample interval of every model parameter.
  ./(2.^paragenelen - 1);
startpointer = cumsum(paragenelen) - paragenelen;                              % the position pointer of the starting point when do model parameter loop.
totalgenelen = sum(paragenelen);                                               % the total length of the binary gene codes of all model parameters.

modelpopuprevcode = NaN*ones(popunum, totalgenelen);                           % the binary codes of all the previous model population, every row respresents a sequence of binary code of a model.
fitnessprev = zeros(1, popunum);                                               % the fitness of every model from modelpopuprevcode.

modelpopunewcode = round(rand(popunum, totalgenelen));                         % the binary codes of all the new model population, every row respresents a sequence of binary code of a model from modelpopunew.
modelpopunew = NaN*ones(paranum, popunum);                                     % the new model population of the 1st or current generation, every column respresents a model.
for i = 1:1:popunum
  modelpopunew(:, i) = DeCodeGA(pararange(:, 1), parasmplintv, paragenelen, ...
    modelpopunewcode(i, :));
end

ign = 0;                                                                       % the ignth generation model population.

while(1)
  ign = ign + 1;
  Tempupdate = Tempupdate*0.95;
  Tempselect = Tempselect*0.95;
  % Fitness evaluation
  fitnessnew = EvaluateFitness(fwdfun, modelpopunew, obsvstn, obsvdat, ...     % the fitness of every model from modelpopunew.
    fitnessmaxlimit);
  % Stopping criteria
  [fitnessmax, ifitmax] = max(fitnessnew);                                     % fitnessmax is the maximum value of fitnessinit, ifitmax is the index of fitnessmax.
  fitnessave = mean(fitnessnew);                                               % the average value of fitnessinit.
  if(fitnessmax == 0)
    error(['The function StandardGeneticAlgorith goes wrong, because all ', ...
      'model population have a very bad fitness. Please try again.']);
  end
  if(fitnessmax > allowquitfitmax ...
    || fitnessave/fitnessmax > allowquitfitavemaxratio)
    fprintf(['\nThe generation number of genetic heredity iteration is: ', ...
      '%d.\n'], ign);
    break;
  end
  if(ign >= ignmax)
    warning('MATLAB:MICMaximumIterationTime', ['The generation number of ', ...
      'genetic heretidy iteration comes to %d, and force the Genetic ', ...
      'Algorithm loop to exit. So the solution may be not enough ', ...
      'accurate.'], ign);
    break;
  end
  % Updating population for a generation
  for i = 1:1:popunum
    if(fitnessprev(i) > fitnessnew(i))
      probupdate = exp((fitnessnew(i) - fitnessprev(i))/Tempupdate);           % the probability of a old model from modelpopuprev being not updated, i.e. being used as a new model in modelpopunew.
      if(rand < probupdate)
        modelpopunewcode(i, :) = modelpopuprevcode(i, :);
        fitnessnew(i) = fitnessprev(i);
      end
    end
  end
  modelpopuprevcode = modelpopunewcode;
  fitnessprev = fitnessnew;
  % Selection (Fitness proportionate selection)
  probselect = exp(fitnessprev./Tempselect) ...                                % the probability of a model from modelpopuinit being selected to be a father of mother.
    ./sum(exp(fitnessprev./Tempselect));
  cumprobselect = cumsum(probselect);                                          % the cumulative probability density distribution of every model from modelpopuinit.
  for i = 1:1:popunum
    iprob = find(rand <= cumprobselect, 1, 'first');                           % the iprobth probability interval is selected, i.e. the iprobth row of modelpopuinitcode is selected to be a new model binary code sequence.
    modelpopunewcode(i, :) = modelpopuprevcode(iprob, :);
  end
  % Crossover (Multi-point)
  for i = 1:2:popunum
    if(rand <= probcrossover)
      for j = 1:paranum
        crossoverpoint = randi([1, paragenelen(j)]);                           % the crossover point position of the jth model parameter.
        crossoverlowerbound = startpointer(j) + crossoverpoint + 1;            % the lower bound of the gene fragment to crossover.
        crossoverupperbound = startpointer(j) + paragenelen(j);                % the upper bound of the gene fragment to crossover.
        tempgenefragment(1:paragenelen(j) - crossoverpoint) ...                % the temporary gene fragment from the ith column of modelpopunewcode.
          = modelpopunewcode(i, crossoverlowerbound:crossoverupperbound);
        modelpopunewcode(i, crossoverlowerbound:crossoverupperbound) ...
          = modelpopunewcode(i + 1, crossoverlowerbound:crossoverupperbound);
        modelpopunewcode(i + 1, crossoverlowerbound:crossoverupperbound) ...
          = tempgenefragment(1:paragenelen(j) - crossoverpoint);
      end
    end
  end
  % Mutation
  for i = 1:1:popunum
    if(rand <= probmutation)
      mutationpoint = randi([1, totalgenelen]);                                % the mutation point position of the ith row model binary code sequence of modelpopunewcode.
      modelpopunewcode(i, mutationpoint) = 1 ...
        - modelpopunewcode(i, mutationpoint);
    end
  end
  % Decoding
  for i = 1:1:popunum
    modelpopunew(:, i) = DeCodeGA(pararange(:, 1), parasmplintv, ...
      paragenelen, modelpopunewcode(i, :));
  end
end

msltn = modelpopunew(:, ifitmax);
rsdl2norm = norm(fwdfun(msltn, obsvstn) - obsvdat, 2);

end

function modelfitness = EvaluateFitness(fwdfun, modelset, obsvstn, obsvdat, ...
  fitnessmaxlimit)
%
% modelfitness = EvaluateFitness(fwdfun, modelset, obsvstn, obsvdat, ...
%   fitnessmaxlimit)
%
% This is a program used to evaluate model fitness by a genetic
% algorithm program.
% Written by Tche L. at USTC, 2015/12.

% modelfitness: a vector whose size is [1, modelnum],
%               the fitnesses of these models from modelset.

% fwdfun: a function handle,
%         the forward function of the inversion problem of
%         the genetic algorithm.
% modelset: a matrix whose size is [paranum, modelnum],
%           the model set whose fitness we want to evaluate.
% obsvstn: a vector whose size is [datanum, ~],
%          the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum, 1],
%          the observe data used to inverse the problem
%          obsvdat = fwdfun(model,obsvstn).
% fitnessmaxlimit: a constant variable,
%                  the maximum value of fitness that can be tolerated.

[~, modelnum] = size(modelset);                                                % the model number of modelset.

modelfitness = zeros(1, modelnum);

for i = 1:1:modelnum
  prdcdat = fwdfun(modelset(:, i), obsvstn);                                   % the predict data when using the ith column of modelset as a model.
  rsdl2norm = norm(obsvdat - prdcdat, 2);
  if(rsdl2norm < fitnessmaxlimit)
    modelfitness(i) = (fitnessmaxlimit - rsdl2norm)/fitnessmaxlimit;
  else
    modelfitness(i) = 0;
  end
end

end

function model = DeCodeGA(paramin, parasmplintv, paragenelen, modelcode)
%
% model = DeCodeGA(paramin, parasmplintv, paragenelen, modelcode)
%
% This is a binary decoding program used by a Genetic Algorithm program.
%
% Written by Tche L. at USTC, 2015/12.

% model: a vector whose size si [paranum, 1],
%        the decoding result model.

% paramin: a vector whose size is [paranum, 1],
%          the available minimum value of every model parameter.
% parasmplintv: a vector whose size is [paranum, 1],
%               the sample interval of every model parameter.
% paragenelen: a vector whose size is [paranum, 1],
%              the length of the gene fragment of every model parameter.
% modelcode: a vector whose size is [1, totalgenelen],
%            the model binary gene code sequence we want to decode.

paranum = length(paragenelen);                                                 % the number of model parameters.

startpointer = 0;                                                              % the position pointer of the starting point when do model parameter loop.
paramultiple = zeros(paranum, 1);                                              % the multiple number of the difference between model and paramin to parasmplintv.
for i = 1:1:paranum
  for j = 1:1:paragenelen(i)
    num2exp = 2^(paragenelen(i) - j);                                          % the exponent of 2 corresponding to the jth binary bit from high to low bit.
    paramultiple(i) = paramultiple(i) + modelcode(startpointer + j)*num2exp;
  end
  startpointer = startpointer + paragenelen(i);
end
model = paramin + paramultiple.*parasmplintv;

end
