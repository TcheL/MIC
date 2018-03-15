function [msltn,rsdl2norm] = StandardGeneticAlgorithm(fwdfun,pararange,paragenelen,obsvstn,obsvdat,fitnessmaxlimit,allowquitfitmax,allowquitfitavemaxratio)
% [msltn,rsdl2norm] = StandardGeneticAlgorithm(fwdfun,pararange,paragenelen,obsvstn,obsvdat,fitnessmaxlimit,allowquitfitmax,allowquitfitavemaxratio)
% This is a standard genetic algorithm inversion program for the problem obsvdat = fwdfun(msltn,obsvstn) in the range restrain pararange of model parameters.
% Written and recommended by Tche.L. from USTC, 2015.12.
%
% msltn: a vector whose size is [paranum,1], the model solution given by this program.
% rsdl2norm: a constant variable, the 2-norm of the residual for the model msltn.
%
% fwdfun: a function handle, the forward function.
% pararange: a matrix whose size is [paranum,2], every row is corresponding to a model parameter, the 1st column is its minimum and the 2nd column is its maximum.
% paragenelen: a vector whose size is [paranum,1], the length of the gene fragment of every model parameter.
% obsvstn: a vector whose size is [datanum,~], the observe station point x for the observe data.
% obsvdat: a vector whose size is [datanum,1], the observe data used to inverse the problem obsvdat = fwdfun(msltn,obsvstn).
% fitnessmaxlimit: a constant variable, the maximum value of fitness that can be tolerated.
% allowquitfitmax: a constant variable who is in [0,1], the fitness value when the maximum fitness of a generation population is greater than it we can quit the genetic heredity iteration.
% allowquitfitavemaxratio: a constant variable who is in [0,1], the ratio value when the ratio of the average fitness to the maximum fitness of a generation population is greater than it we can quit the genetic heredity iteration.

paranum = length(paragenelen);                                                  % the number of model parameters.
popunum = 2^10;                                                                 % the population number of every generation.
probcrossover = 0.7;                                                            % the probability of crossover.
probmutation = 0.05;                                                            % the probability of mutation.
ignmax = 100;                                                                   % the generation number of genetic heredity.

% Encoding (Coding of model parameters)
parasmplintv = (pararange(:,2) - pararange(:,1)) ...                            % the sample interval of every model parameter.
    ./(2.^paragenelen - 1);
startpointer = cumsum(paragenelen) - paragenelen;                           	% the position pointer of the starting point when do model parameter loop.
totalgenelen = sum(paragenelen);                                                % the total length of the binary gene codes of all model parameters.
modelpopuinitcode = round(rand(popunum,totalgenelen));                          % the binary codes of all the initial model population, every row respresents a sequence of binary code of a model from modelpopuinit.
modelpopuinit = NaN*ones(paranum,popunum);                                      % the initial model population of the 1st or current generation, every column respresents a model.
for i = 1:1:popunum
    modelpopuinit(:,i) = DeCodeGA(pararange(:,1),parasmplintv,paragenelen,modelpopuinitcode(i,:));
end

ign = 0;                                                                        % the ignth generation model population.
r2normset = Inf*ones(1,popunum);                                                % the set of the 2-norm of the residual for all models of the modelpopuinit.
modelfitness = zeros(1,popunum);                                                % the fitness of every model from modelpopuinit.
modelpopunew = NaN*ones(paranum,popunum);                                       % the new model population of the next generation, every column respresents a model.
modelpopunewcode = NaN*ones(popunum,totalgenelen);                              % the binary codes of all the new model population, every row respresents a sequence of binary code of a model from modelpopunew.
while(1)
    ign = ign + 1;
    % Fitness evaluation
    for i = 1:1:popunum
        prdcdat = fwdfun(modelpopuinit(:,i),obsvstn);                           % the predict data when using the ith column of modelpopuinit as a model.
        r2normset(i) = norm(obsvdat - prdcdat,2);
        if(r2normset(i) < fitnessmaxlimit)
            modelfitness(i) = (fitnessmaxlimit - r2normset(i))/fitnessmaxlimit;
        else
            modelfitness(i) = 0;
        end
    end
    % Stopping criteria
    [fitnessmax,ifitmax] = max(modelfitness);                                   % fitnessmax is the maximum value of modelfitness, ifitmax is the index of fitnessmax.
    fitnessave = mean(modelfitness);                                            % the average value of modelfitness.
    if(fitnessmax == 0)
        error('The function StandardGeneticAlgorith goes wrong, because all model population have a very bad fitness. Please try again.');
    end
    if(fitnessmax > allowquitfitmax || fitnessave/fitnessmax > allowquitfitavemaxratio)
        fprintf('\nThe generation number of genetic heredity iteration is: %d.\n',ign);
        break;
    end
    if(ign >= ignmax)
        warning('MATLAB:MyLargestIterationTimes','The generation number of genetic heretidy iteration comes to %d, and force the Genetic Algorithm loop to exit. So the solution may be not enough accurate.',ign);
        break;
    end
    % Selection (Fitness proportionate selection)
    probselect = modelfitness./sum(modelfitness);                               % the probability of a model from modelpopuinit being selected to be a father of mother.
    cumprobselect = cumsum(probselect);                                         % the cumulative probability density distribution of every model from modelpopuinit.
    for i = 1:1:popunum
        iprob = find(rand <= cumprobselect,1,'first');                          % the iprobth probability interval is selected, i.e. the iprobth row of modelpopuinitcode is selected to be a new model binary code sequence.
        modelpopunewcode(i,:) = modelpopuinitcode(iprob,:);
    end
    % Crossover (Multi-point)
    for i = 1:2:popunum
        if(rand <= probcrossover)
            for j = 1:paranum
                crossoverpoint = randi([1,paragenelen(j)]);                   	% the crossover point position of the jth model parameter.
                crossoverlowerbound = startpointer(j) + crossoverpoint + 1;  	% the lower bound of the gene fragment to crossover.
                crossoverupperbound = startpointer(j) + paragenelen(j);       	% the upper bound of the gene fragment to crossover.
                tempgenefragment(1:paragenelen(j) - crossoverpoint) = ...    	% the temporary gene fragment from the ith column of modelpopunewcode.
                    modelpopunewcode(i,crossoverlowerbound:crossoverupperbound);
                modelpopunewcode(i,crossoverlowerbound:crossoverupperbound) = modelpopunewcode(i + 1,crossoverlowerbound:crossoverupperbound);
                modelpopunewcode(i + 1,crossoverlowerbound:crossoverupperbound) = tempgenefragment(1:paragenelen(j) - crossoverpoint);
            end
        end
    end
    % Mutation
    for i = 1:1:popunum
        if(rand <= probmutation)
            mutationpoint = randi([1,totalgenelen]);                            % the mutation point position of the ith row model binary code sequence of modelpopunewcode.
            modelpopunewcode(i,mutationpoint) = 1 - modelpopunewcode(i,mutationpoint);
        end
    end
    % Decoding
    for i = 1:1:popunum
        modelpopunew(:,i) = DeCodeGA(pararange(:,1),parasmplintv,paragenelen,modelpopunewcode(i,:));
    end
    % Iterating to do loop
    modelpopuinit = modelpopunew;
    modelpopuinitcode = modelpopunewcode;
end

msltn = modelpopuinit(:,ifitmax);
rsdl2norm = r2normset(ifitmax);

end

function model = DeCodeGA(paramin,parasmplintv,paragenelen,modelcode)
% model = DeCodeGA(paramin,parasmplintv,paragenelen,modelcode)
% This is a binary decoding program used by a Genetic Algorithm program.
% Written by Tche.L. from USTC, 2015,12.

% model: a vector whose size si [paranum,1], the decoding result model.

% paramin: a vector whose size is [paranum,1], the available minimum value of every model parameter.
% parasmplintv: a vector whose size is [paranum,1], the sample interval of every model parameter.
% paragenelen: a vector whose size is [paranum,1], the length of the gene fragment of every model parameter.
% modelcode: a vector whose size is [1,totalgenelen], the model binary gene code sequence we want to decode.

paranum = length(paragenelen);                                                  % the number of model parameters.

startpointer = 0;                                                               % the position pointer of the starting point when do model parameter loop.
paramultiple = zeros(paranum,1);                                                % the multiple number of the difference between model and paramin to parasmplintv.
for i = 1:1:paranum
    for j = 1:1:paragenelen(i)
        num2exp = 2^(paragenelen(i) - j);                                       % the exponent of 2 corresponding to the jth binary bit from high to low bit.
        paramultiple(i) = paramultiple(i) + modelcode(startpointer + j)*num2exp;
    end
    startpointer = startpointer + paragenelen(i);
end
model = paramin + paramultiple.*parasmplintv;

end