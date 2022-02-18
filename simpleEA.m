function [bestSoFarFit ,bestSoFarSolution ...
    ]=simpleEA( ...  % name of your simple EA function
    fitFunc, ... % name of objective/fitness function
    T, ... % total number of evaluations
    input) % replace it by your input arguments

% Check the inputs
if isempty(fitFunc)
  warning(['Objective function not specified, ''' objFunc ''' used']);
  fitFunc = 'objFunc';
end
if ~ischar(fitFunc)
  error('Argument FITFUNC must be a string');
end
if isempty(T)
  warning(['Budget not specified. 1000000 used']);
  T = '1000000';
end
eval(sprintf('objective=@%s;',fitFunc));
% Initialise variables
nbGen = 0; % generation counter
nbEval = 0; % evaluation counter
bestSoFarFit = 0; % best-so-far fitness value
bestSoFarSolution = NaN; % best-so-far solution
p_size=4
%recorders
fitness_gen=[0]; % record the best fitness so far
solution_gen=[];% record the best phenotype of each generation
fitness_pop=[];% record the best fitness in current population 
%% Below starting your code

% Initialise a population
%% TODO
P=fix(rand(1,p_size)*31);
B=dec2bin(P);
% Evaluate the initial population
%% TODO

% Start the loop
while (nbEval<T) 
% Reproduction (selection, crossver)
%% TODO
%selection
%eval
evals=objective(P);
nbEval=nbEval+p_size
prob=evals/sum(evals);
probs=prob
len_gen=5
for i=2:p_size
    probs(i)=prob(i)+probs(i-1);
end
parents=[];
for i=1:2
    p=rand(1);
    probs_aux=[probs,p];
    idx=find(sort(probs_aux)==p);
    parents=[parents P(idx)];
end
parents_B=dec2bin(parents);

%crossover
crossover_idx=fix(rand(1,1)*len_gen)+1;
offs1=[parents_B(1,1:crossover_idx),parents_B(2,crossover_idx+1:len_gen)]
offs2=[parents_B(2,1:crossover_idx),parents_B(1,crossover_idx+1:len_gen)]

% Mutation
mu_idx=fix(rand(1,1)*len_gen)+1;
offs1(mu_idx)=num2str(1-str2num(offs1(mu_idx)));
mu_idx=fix(rand(1,1)*len_gen)+1;
offs2(mu_idx)=num2str(1-str2num(offs2(mu_idx)));
offs1=bin2dec(offs1);
offs2=bin2dec(offs2);
P_next=[P offs1 offs2];
evals=objective(P_next);
P_next=[P_next;evals]'
P_next=sortrows(P_next,-2) 
P=P_next(1:p_size,1)

x=1

%% TODO
solution_gen=[solution_gen P_next(1,1)];% record the best phenotype of each generation
fitness_pop=[fitness_pop P_next(1,2)]
if P_next(1,2)>fitness_gen(length(fitness_gen))
    bestSoFarFit=[fitness_gen P_next(1,2)]
else
    fitness_gen=[fitness_gen fitness_gen(fitness_gen(length(fitness_gen)))]
end
nbEval=nbEval+p_size+2
bestSoFarFit=
end





