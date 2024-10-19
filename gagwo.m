
nTasks = 10;                % Number of Tasks
nMachines = 5;             % Number of Machines
VarSize = [1 nTasks];      % Size of  task assigned to a machine
VarMin = 1;                % Lower Bound of Variables 
VarMax = nMachines;        % Upper Bound of Variables 

% GAGWO Parameters
MaxIt = 100;              % Maximum Number of Iterations
nPop = 30;                % Population Size
Pc = 0.7;                 % Crossover Probability
Pm = 0.3;                 % Mutation Probability
beta = 1.5;               % Mutation Rate
nCrossover = round(Pc * nPop / 2) * 2;  % Number of Crossover Pairs
nMutation = round(Pm * nPop);           % Number of Mutants

% Grey Wolf Optimizer Parameters
alpha = struct('Position', [], 'Cost', inf);
betaStruct = struct('Position', [], 'Cost', inf);
delta = struct('Position', [], 'Cost', inf);

empty_individual.Position = [];
empty_individual.Cost = [];
pop = repmat(empty_individual, nPop, 1);

% Initialize Population
for i = 1:nPop
    pop(i).Position = randi([VarMin, VarMax], VarSize);
    pop(i).Cost = ObjectiveFunction(pop(i).Position);
    % Update Alpha, betaStruct, Delta
    if pop(i).Cost < alpha.Cost
        delta = betaStruct;
        betaStruct = alpha;
        alpha = pop(i);
    elseif pop(i).Cost < betaStruct.Cost
        delta = betaStruct;
        betaStruct = pop(i);
    elseif pop(i).Cost < delta.Cost
        delta = pop(i);
    end
end

% GAGWO Main Loop
BestCosts = nan(MaxIt, 1);

for it = 1:MaxIt
    
   
    a = 2 - it * (2 / MaxIt);
    
    for i = 1:nPop
        A1 = 2 * a * rand(VarSize) - a;
        C1 = 2 * rand(VarSize);
        D_alpha = abs(C1 .* alpha.Position - pop(i).Position);
        X1 = alpha.Position - A1 .* D_alpha;
        
        A2 = 2 * a * rand(VarSize) - a;
        C2 = 2 * rand(VarSize);
        D_beta = abs(C2 .* betaStruct.Position - pop(i).Position);
        X2 = betaStruct.Position - A2 .* D_beta;
        
        A3 = 2 * a * rand(VarSize) - a;
        C3 = 2 * rand(VarSize);
        D_delta = abs(C3 .* delta.Position - pop(i).Position);
        X3 = delta.Position - A3 .* D_delta;
        
        % Update Position (Discrete Assignment)
        newPosition = round((X1 + X2 + X3) / 3);
        newPosition = max(newPosition, VarMin);
        newPosition = min(newPosition, VarMax);
        pop(i).Position = newPosition;
        
        
        pop(i).Cost = ObjectiveFunction(pop(i).Position);
        
        % Update Alpha, betaStruct, Delta
        if pop(i).Cost < alpha.Cost
            delta = betaStruct;
            betaStruct = alpha;
            alpha = pop(i);
        elseif pop(i).Cost < betaStruct.Cost
            delta = betaStruct;
            betaStruct = pop(i);
        elseif pop(i).Cost < delta.Cost
            delta = pop(i);
        end
    end
    
    % Crossover
    popc = repmat(empty_individual, nCrossover / 2, 2);
    for k = 1:nCrossover / 2
        i1 = randi([1 nPop]);
        i2 = randi([1 nPop]);
        
        p1 = pop(i1);
        p2 = pop(i2);
        
        [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position);
        
        popc(k, 1).Cost = ObjectiveFunction(popc(k, 1).Position);
        popc(k, 2).Cost = ObjectiveFunction(popc(k, 2).Position);
    end
    popc = popc(:);
    
    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation
        i = randi([1 nPop]);
        p = pop(i);
        
        popm(k).Position = Mutate(p.Position, beta, VarMin, VarMax);
        popm(k).Cost = ObjectiveFunction(popm(k).Position);
    end
    
    % Merge
    pop = [pop
           popc
           popm];
       
    % Sort Population
    [~, SortOrder] = sort([pop.Cost]);
    pop = pop(SortOrder);
    
    % Truncate Extra Members
    pop = pop(1:nPop);
    
    % Store Best Cost
    BestCosts(it) = alpha.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    
end


figure;
plot(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


function cost = ObjectiveFunction(schedule)
    
    nMachines = 5; % Example number of machines
    taskTimes = [2, 4, 6, 3, 7, 5, 8, 4, 9, 3]; % Example task processing times
    machineTimes = zeros(1, nMachines);
    
    for i = 1:length(schedule)
        machine = schedule(i);
        machineTimes(machine) = machineTimes(machine) + taskTimes(i);
    end
    
    
    cost = max(machineTimes);
end

% Crossover Function
function [y1, y2] = Crossover(x1, x2)
    alpha = rand(size(x1));
    y1 = round(alpha.*x1 + (1-alpha).*x2);
    y2 = round(alpha.*x2 + (1-alpha).*x1);
end


function y = Mutate(x, mu, VarMin, VarMax)
    nVar = numel(x);
    nMu = ceil(mu * nVar);
    
    j = randperm(nVar);
    nMu = min(nMu, nVar);
j = j(1:nMu);
    
    y = x;
    y(j) = randi([VarMin, VarMax], size(j));
    
    y = max(y, VarMin);
    y = min(y, VarMax);
end
