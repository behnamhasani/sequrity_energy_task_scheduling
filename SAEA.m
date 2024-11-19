function [makespan, loadBalance, avgSecurityLevel, totalEnergyCost] = SAEA(tasks, VMs, maxIterations)
    
    numTasks = length(tasks.executionTime);
    numVMs = length(VMs.capacity.CPU);
    numSquirrels = 50; 
    
    % Initialize squirrel population
    squirrels = initializeSquirrels(numTasks, numVMs, numSquirrels);
    
    % Initialize metrics
    bestMakespan = inf;
    bestLoadBalance = 0;
    totalSecurityMetric = 0; 
    totalEnergyCost = 0;

    % Main loop for optimization
    for iter = 1:maxIterations
       
        fitness = evaluateFitness_SAEA(squirrels, tasks, VMs);
        
        % Perform parallel search and migration
        squirrels = parallelSearch(squirrels, fitness, numTasks, numVMs);
        
        % Communication step to share best solutions
        squirrels = communicateBestSolutions(squirrels);
        
        % Find the best squirrel in this iteration
        [~, bestIdx] = min(fitness);
        bestSquirrel = squirrels(bestIdx, :);
        
        currentMakespan = CalculateMakespan(bestSquirrel, tasks, VMs);
        currentLoadBalance = CalculateLoadBalance(bestSquirrel, tasks, VMs);
        currentEnergy = CalculateEnergyConsumption(bestSquirrel, tasks, VMs);
        currentSecurity = CalculateSecurityBias(bestSquirrel, tasks, VMs);
        
        % Update overall metrics
        bestMakespan = min(bestMakespan, currentMakespan);
        totalEnergyCost = totalEnergyCost + currentEnergy;
        totalSecurityMetric = totalSecurityMetric + currentSecurity;
        bestLoadBalance = max(bestLoadBalance, currentLoadBalance);
        
        fprintf('Iteration %d: Best Fitness = %.2f, Makespan = %.2f\n', iter, fitness(bestIdx), currentMakespan);
    end
    
    makespan = bestMakespan;
    loadBalance = bestLoadBalance;
    avgSecurityLevel = totalSecurityMetric / maxIterations;
    totalEnergyCost = totalEnergyCost;
end



function squirrels = initializeSquirrels(numTasks, numVMs, numSquirrels)
    squirrels = randi([1, numVMs], numSquirrels, numTasks);  
end

% Evaluate Fitness 
function fitness = evaluateFitness_SAEA(squirrels, tasks, VMs)
    numSquirrels = size(squirrels, 1);
    fitness = zeros(numSquirrels, 1);
    
    for i = 1:numSquirrels
        taskAssignment = squirrels(i, :);
        makespan = CalculateMakespan(taskAssignment, tasks, VMs);
        energy = CalculateEnergyConsumption(taskAssignment, tasks, VMs);
        security = CalculateSecurityBias(taskAssignment, tasks, VMs);
        loadBalance = CalculateLoadBalance(taskAssignment, tasks, VMs);
        
        % Combine objectives 
        weights = struct('makespan', 0.4, 'energy', 0.3, 'security', 0.2, 'loadBalance', 0.1);
        fitness(i) = weights.makespan * makespan + ...
                     weights.energy * energy + ...
                     weights.security * security + ...
                     weights.loadBalance * loadBalance;
    end
end

function makespan = CalculateMakespan(taskAssignment, tasks, VMs)
    vmCompletionTimes = zeros(1, length(VMs.capacity.CPU));
    for task = 1:length(taskAssignment)
        vm = taskAssignment(task);
        adjustedTime = tasks.executionTime(task) * (tasks.resource.CPU(task) / VMs.capacity.CPU(vm));
        vmCompletionTimes(vm) = vmCompletionTimes(vm) + adjustedTime;
    end
    makespan = max(vmCompletionTimes);
end

function energy = CalculateEnergyConsumption(taskAssignment, tasks, VMs)
    energy = 0;
    for vm = 1:length(VMs.capacity.CPU)
        assignedTasks = find(taskAssignment == vm);
        if ~isempty(assignedTasks)
            T_busy = sum(tasks.executionTime(assignedTasks));
            energy_busy = VMs.power.busy(vm) * T_busy;
            E_security = sum(tasks.sensitivity.confidentiality(assignedTasks)) + ...
                         sum(tasks.sensitivity.integrity(assignedTasks)) + ...
                         sum(tasks.sensitivity.authentication(assignedTasks));
            energy = energy + energy_busy + E_security;
        else
            T_idle = 200; 
            energy_idle = VMs.power.idle(vm) * T_idle;
            energy = energy + energy_idle;
        end
    end
end

function security = CalculateSecurityBias(taskAssignment, tasks, VMs)
    security = 0;
    for task = 1:length(taskAssignment)
        vm = taskAssignment(task);
        confidentialityOverhead = tasks.sensitivity.confidentiality(task) * (1 - VMs.trustLevel(vm));
        integrityOverhead = tasks.sensitivity.integrity(task) * (1 - VMs.trustLevel(vm));
        authenticationOverhead = tasks.sensitivity.authentication(task) * (1 - VMs.trustLevel(vm));
        security = security + (confidentialityOverhead + integrityOverhead + authenticationOverhead);
    end
end

function loadBalance = CalculateLoadBalance(taskAssignment, tasks, VMs)
    vmLoad.CPU = zeros(1, length(VMs.capacity.CPU));
    for vm = 1:length(VMs.capacity.CPU)
        assignedTasks = find(taskAssignment == vm);
        vmLoad.CPU(vm) = sum(tasks.resource.CPU(assignedTasks));
    end
    avgCpuUtilization = mean(vmLoad.CPU ./ VMs.capacity.CPU);
    cpuImbalance = sum((vmLoad.CPU ./ VMs.capacity.CPU - avgCpuUtilization).^2);
    loadBalance = 1 / (1 + cpuImbalance);
end

function squirrels = parallelSearch(squirrels, fitness, numTasks, numVMs)
    numSquirrels = size(squirrels, 1);
    [~, sortedIdx] = sort(fitness); 
    group1 = squirrels(sortedIdx(1:round(numSquirrels / 2)), :);  % Best group
    group2 = squirrels(sortedIdx(round(numSquirrels / 2)+1:end), :);  % Worst group
    
    % Migrate worst group towards best group
    for i = 1:size(group2, 1)
        leader = group1(randi([1, size(group1, 1)]), :);  
        group2(i, :) = migrate(group2(i, :), leader, numTasks, numVMs);
    end
    
    squirrels = [group1; group2];  
end

% Migration Function
function migratedSquirrel = migrate(squirrel, leader, numTasks, numVMs)
    migrationRate = 0.5;  
    
    for i = 1:numTasks
        if rand() < migrationRate
            squirrel(i) = leader(i);  % Move task assignment towards leader
        end
    end
    
    migratedSquirrel = squirrel;
end

% Communication Step
function squirrels = communicateBestSolutions(squirrels)
    numSquirrels = size(squirrels, 1);
    bestSquirrel = squirrels(1, :);  
    
    for i = 1:numSquirrels
        if rand() < 0.1  
            squirrels(i, :) = bestSquirrel;
        end
    end
end
