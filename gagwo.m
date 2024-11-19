function [makespan, loadBalance, avgSecurityLevel, totalEnergyCost] = gagwo(tasks, VMs, maxIterations)
    numTasks = length(tasks.executionTime);
    numVMs = length(VMs.capacity.CPU);
    VarSize = numTasks;      
    VarMin = 1;              % Minimum Value (VM Index)
    VarMax = numVMs;         % Maximum Value (VM Index)
    nPop = 30;               % Population Size
    Pc = 0.7;                % Crossover Probability
    Pm = 0.3;                % Mutation Probability
    beta = 1.5;              % Mutation Rate
    nCrossover = round(Pc * nPop / 2) * 2;  % Number of Crossover Pairs
    nMutation = round(Pm * nPop);           % Number of Mutants

    empty_individual.Position = [];
    empty_individual.Cost = inf;
    empty_individual.Makespan = inf;
    empty_individual.Energy = inf;
    empty_individual.Security = inf;
    empty_individual.LoadBalance = inf;
    pop = repmat(empty_individual, nPop, 1);

    % Initialize Metrics
    bestMakespan = inf;
    totalEnergyCost = 0;
    totalSecurityMetric = 0;
    bestLoadBalance = 0;

    % Initialize Population
    for i = 1:nPop
        pop(i).Position = randi([VarMin, VarMax], [1, VarSize]);  
        [pop(i).Cost, pop(i).Makespan, pop(i).Energy, pop(i).Security, pop(i).LoadBalance] = ...
            EvaluateFitness(pop(i).Position, tasks, VMs);  
    end

    % Main GAGWO Loop
    for iter = 1:maxIterations
        for i = 1:nPop
            [pop(i).Cost, makespan, energy, security, loadBalance] = EvaluateFitness(pop(i).Position, tasks, VMs);

            bestMakespan = min(bestMakespan, makespan);
            totalEnergyCost = totalEnergyCost + energy;
            totalSecurityMetric = totalSecurityMetric + security;
            bestLoadBalance = max(bestLoadBalance, loadBalance);
        end

        popCrossover = repmat(empty_individual, nCrossover, 1);
        for k = 1:2:nCrossover
            i1 = randi(nPop);
            i2 = randi(nPop);

            [popCrossover(k).Position, popCrossover(k + 1).Position] = ...
                Crossover(pop(i1).Position, pop(i2).Position);

            [popCrossover(k).Cost, popCrossover(k).Makespan, popCrossover(k).Energy, ...
                popCrossover(k).Security, popCrossover(k).LoadBalance] = ...
                EvaluateFitness(popCrossover(k).Position, tasks, VMs);

            [popCrossover(k + 1).Cost, popCrossover(k + 1).Makespan, ...
                popCrossover(k + 1).Energy, popCrossover(k + 1).Security, ...
                popCrossover(k + 1).LoadBalance] = ...
                EvaluateFitness(popCrossover(k + 1).Position, tasks, VMs);
        end

        popMutation = repmat(empty_individual, nMutation, 1);
        for k = 1:nMutation
            i = randi(nPop);
            popMutation(k).Position = Mutate(pop(i).Position, beta, VarMin, VarMax);
            [popMutation(k).Cost, popMutation(k).Makespan, popMutation(k).Energy, ...
                popMutation(k).Security, popMutation(k).LoadBalance] = ...
                EvaluateFitness(popMutation(k).Position, tasks, VMs);
        end

        pop = [pop; popCrossover; popMutation];

        [~, sortIdx] = sort([pop.Cost]);
        pop = pop(sortIdx);

        % Retain Top Individuals
        pop = pop(1:nPop);
    end

    makespan = bestMakespan;
    loadBalance = bestLoadBalance;
    avgSecurityLevel = totalSecurityMetric / maxIterations;
    totalEnergyCost = totalEnergyCost;
end

function [fitness, makespan, energy, security, loadBalance] = EvaluateFitness(positionMatrix, tasks, VMs)
    makespan = CalculateMakespan(positionMatrix, tasks, VMs);
    energy = CalculateEnergyConsumption(positionMatrix, tasks, VMs);
    security = CalculateSecurityBias(positionMatrix, tasks, VMs);
    loadBalance = CalculateLoadBalance(positionMatrix, tasks, VMs);

    wMakespan = 0.4;  
    wEnergy = 0.3;    
    wSecurity = 0.2;  
    wLoadBalance = 0.1;  

    fitness = wMakespan * makespan + wEnergy * energy + wSecurity * security - wLoadBalance * loadBalance;
end

function makespan = CalculateMakespan(positionMatrix, tasks, VMs)
    vmCompletionTimes = zeros(1, length(VMs.capacity.CPU));  
    for vm = 1:length(VMs.capacity.CPU)
        assignedTasks = find(positionMatrix == vm); 
        for task = assignedTasks
            adjustedExecutionTime = tasks.executionTime(task) * (tasks.resource.CPU(task) / VMs.capacity.CPU(vm));
            vmCompletionTimes(vm) = vmCompletionTimes(vm) + adjustedExecutionTime;
        end
    end
    makespan = max(vmCompletionTimes); 
end

function energy = CalculateEnergyConsumption(positionMatrix, tasks, VMs)
    energy = 0;
    for vm = 1:length(VMs.capacity.CPU)
        assignedTasks = find(positionMatrix == vm);
        T_busy = sum(tasks.executionTime(assignedTasks));
        energy_busy = VMs.power.busy(vm) * T_busy;
        energy = energy + energy_busy;
    end
end

function security = CalculateSecurityBias(positionMatrix, tasks, VMs)
    security = 0;
    for task = 1:length(tasks.executionTime)
        vm = positionMatrix(task);
        confidentialityOverhead = tasks.sensitivity.confidentiality(task) * (1 - VMs.trustLevel(vm));
        integrityOverhead = tasks.sensitivity.integrity(task) * (1 - VMs.trustLevel(vm));
        authenticationOverhead = tasks.sensitivity.authentication(task) * (1 - VMs.trustLevel(vm));
        security = security + (confidentialityOverhead + integrityOverhead + authenticationOverhead);
    end
end

function loadBalance = CalculateLoadBalance(positionMatrix, tasks, VMs)
    vmLoad.CPU = zeros(1, length(VMs.capacity.CPU));
    for vm = 1:length(VMs.capacity.CPU)
        assignedTasks = find(positionMatrix == vm);
        vmLoad.CPU(vm) = sum(tasks.resource.CPU(assignedTasks));
    end
    avgCpuUtilization = mean(vmLoad.CPU ./ VMs.capacity.CPU);
    cpuImbalance = sum((vmLoad.CPU ./ VMs.capacity.CPU - avgCpuUtilization).^2);
    loadBalance = 1 / (1 + cpuImbalance);
end

function [y1, y2] = Crossover(x1, x2)
    alpha = rand(size(x1));
    y1 = round(alpha .* x1 + (1 - alpha) .* x2);
    y2 = round(alpha .* x2 + (1 - alpha) .* x1);
end

function y = Mutate(x, mu, VarMin, VarMax)
    nVar = numel(x);  
    nMu = min(ceil(mu * nVar), nVar);  
    j = randperm(nVar, nMu);  
    y = x;
    y(j) = randi([VarMin, VarMax], size(j));  
end

