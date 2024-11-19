function [bestSolution, bestMakespan, bestEnergy, bestSecurity, bestLoadBalance, iterationResults] = fqaso(tasks, VMs, populationSize, maxIterations)
    
    numTasks = length(tasks.executionTime);
    numVMs = length(VMs.capacity.CPU);
    
    % Initialize population
    positionMatrix = InitializePopulation(numTasks, numVMs, populationSize);
    mappedMatrix = MapMatrix(positionMatrix, numTasks, numVMs);
    
    
    bestSolution = [];
    bestFitness = Inf;
    
    bestMakespan = Inf;
    bestEnergy = Inf;
    bestSecurity = Inf;
    bestLoadBalance = Inf;

    maxObjectives = struct('makespan', -inf, 'energy', -inf, 'security', -inf, 'loadBalance', -inf);
    minObjectives = struct('makespan', inf, 'energy', inf, 'security', inf, 'loadBalance', inf);
    
    iterationResults = struct('iteration', [], 'bestSolution', [], 'bestMakespan', [], 'bestEnergy', [], 'bestSecurity', [], 'bestLoadBalance', []);
    
    % Begin iterations
    for t = 1:maxIterations
        fitnessValues = zeros(1, populationSize);
        for i = 1:populationSize
            % Evaluate fitness
            [fitness, makespan, energy, security, loadBalance] = EvaluateFitness(positionMatrix(:,:,i), tasks, VMs, maxObjectives, minObjectives);
            fitnessValues(i) = fitness;

            % Update the best solution 
            if fitness < bestFitness
                bestFitness = fitness;
                bestSolution = positionMatrix(:,:,i);
                bestMakespan = makespan;
                bestEnergy = energy;
                bestSecurity = security;
                bestLoadBalance = loadBalance;
            end
        end
        
        iterationResults(t).iteration = t;
        iterationResults(t).bestSolution = bestSolution;
        iterationResults(t).bestMakespan = bestMakespan;
        iterationResults(t).bestEnergy = bestEnergy;
        iterationResults(t).bestSecurity = bestSecurity;
        iterationResults(t).bestLoadBalance = bestLoadBalance;
        
        positionMatrix = UpdatePositions(positionMatrix, mappedMatrix, numTasks, numVMs, fitnessValues, bestSolution);
        
        [positionMatrix, mappedMatrix] = MigrateTasks(positionMatrix, mappedMatrix, numTasks, numVMs);
    end
    
    disp('Best solution found:');
    disp(bestSolution);
    disp(['Makespan: ', num2str(bestMakespan)]);
    disp(['Energy Consumption: ', num2str(bestEnergy)]);
    disp(['Security Bias: ', num2str(bestSecurity)]);
    disp(['Load Balance: ', num2str(bestLoadBalance)]);
end


function positionMatrix = InitializePopulation(numTasks, numVMs, populationSize)
    positionMatrix = randi([0, 1], numTasks, numVMs, populationSize); % Binary allocation of tasks to VMs
end

% Map Matrix
function mappedMatrix = MapMatrix(positionMatrix, numTasks, numVMs)
    mappedMatrix = positionMatrix > 0; % Binary map
end

% Evaluate Fitness
function [fitness, makespan, energy, security, loadBalance] = EvaluateFitness(positionMatrix, tasks, VMs, maxObjectives, minObjectives)
    makespan = CalculateMakespan(positionMatrix, tasks, VMs);
    energy = CalculateEnergyConsumption(positionMatrix, tasks, VMs);
    security = CalculateSecurityBias(positionMatrix, tasks, VMs);
    loadBalance = CalculateLoadBalance(positionMatrix, tasks, VMs);

    [maxObjectives, minObjectives] = UpdateMaxMinObjectives(makespan, energy, security, loadBalance, maxObjectives, minObjectives);

    % Normalize objectives
    normalizedMakespan = NormalizeObjective(makespan, minObjectives.makespan, maxObjectives.makespan);
    normalizedEnergy = NormalizeObjective(energy, minObjectives.energy, maxObjectives.energy);
    normalizedLoadBalance = NormalizeObjective(loadBalance, minObjectives.loadBalance, maxObjectives.loadBalance);
    normalizedSecurity = NormalizeObjective(security, minObjectives.security, maxObjectives.security);

    % fuzzy system
    fitness = FuzzySystem(normalizedMakespan, normalizedEnergy, normalizedLoadBalance, normalizedSecurity);
end

%  Makespan
function makespan = CalculateMakespan(positionMatrix, tasks, VMs)
    vmCompletionTimes = zeros(1, length(VMs.trustLevel));  
    for vm = 1:length(VMs.trustLevel)
        assignedTasks = find(positionMatrix(:, vm) == 1); 
        if ~isempty(assignedTasks)
            for task = assignedTasks'
                adjustedExecutionTime = tasks.executionTime(task) * (tasks.resource.CPU(task) / VMs.capacity.CPU(vm));
                vmCompletionTimes(vm) = vmCompletionTimes(vm) + adjustedExecutionTime;
            end
        end
    end
    makespan = max(vmCompletionTimes); 
end

% Energy Consumption
function energy = CalculateEnergyConsumption(positionMatrix, tasks, VMs)
    energy = 0;
    for vm = 1:length(VMs.trustLevel)
        assignedTasks = find(positionMatrix(:, vm) == 1);
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

%  Security Bias
function security = CalculateSecurityBias(positionMatrix, tasks, VMs)
    security = 0;
    for task = 1:length(tasks.executionTime)
        for vm = 1:length(VMs.trustLevel)
            if positionMatrix(task, vm) == 1
                confidentialityOverhead = tasks.sensitivity.confidentiality(task) * (1 - VMs.trustLevel(vm));
                integrityOverhead = tasks.sensitivity.integrity(task) * (1 - VMs.trustLevel(vm));
                authenticationOverhead = tasks.sensitivity.authentication(task) * (1 - VMs.trustLevel(vm));
                security = security + (confidentialityOverhead + integrityOverhead + authenticationOverhead);
            end
        end
    end
end

%  Load Balance
function loadBalance = CalculateLoadBalance(positionMatrix, tasks, VMs)
    vmLoad.CPU = zeros(1, length(VMs.trustLevel));
    for vm = 1:length(VMs.trustLevel)
        assignedTasks = find(positionMatrix(:, vm) == 1); 
        vmLoad.CPU(vm) = sum(tasks.resource.CPU(assignedTasks));  
    end
    avgCpuUtilization = mean(vmLoad.CPU ./ VMs.capacity.CPU);
    cpuImbalance = sum((vmLoad.CPU ./ VMs.capacity.CPU - avgCpuUtilization).^2);
    loadBalance = 1 / (1 + cpuImbalance);
end

% Normalize Objective
function normalizedValue = NormalizeObjective(value, minVal, maxVal)
    if maxVal == minVal
        normalizedValue = 0.5;
    else
        normalizedValue = (value - minVal) / (maxVal - minVal);
    end
    normalizedValue = min(max(normalizedValue, 0), 1);
end

% Fuzzy System
function fitness = FuzzySystem(makespan, energy, loadBalance, security)
    % Create FIS
    fis = mamfis('Name', 'FitnessFIS');
    
    % Input Variables
    fis = addInput(fis, [0 1], 'Name', 'Makespan');
    fis = addInput(fis, [0 1], 'Name', 'Energy');
    fis = addInput(fis, [0 1], 'Name', 'LoadBalance');
    fis = addInput(fis, [0 1], 'Name', 'Security');
    
    %  Output Variable
    fis = addOutput(fis, [0 1], 'Name', 'Fitness');
    
    % Define Membership Functions for Inputs
    fis = addMF(fis, 'Makespan', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Makespan', 'trimf', [0.3 0.5 0.7], 'Name', 'Medium');
    fis = addMF(fis, 'Makespan', 'trapmf', [0.5 0.7 1 1], 'Name', 'High');
    
    fis = addMF(fis, 'Energy', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Energy', 'trimf', [0.3 0.5 0.7], 'Name', 'Medium');
    fis = addMF(fis, 'Energy', 'trapmf', [0.5 0.7 1 1], 'Name', 'High');
    
    fis = addMF(fis, 'LoadBalance', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'LoadBalance', 'trimf', [0.3 0.5 0.7], 'Name', 'Medium');
    fis = addMF(fis, 'LoadBalance', 'trapmf', [0.5 0.7 1 1], 'Name', 'High');
    
    fis = addMF(fis, 'Security', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Security', 'trimf', [0.3 0.5 0.7], 'Name', 'Medium');
    fis = addMF(fis, 'Security', 'trapmf', [0.5 0.7 1 1], 'Name', 'High');
    
    % Define Membership Functions for Output
    fis = addMF(fis, 'Fitness', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Fitness', 'trimf', [0.3 0.5 0.7], 'Name', 'Medium');
    fis = addMF(fis, 'Fitness', 'trapmf', [0.5 0.7 1 1], 'Name', 'High');
    
    % Define Rules
    rules = [
        % Makespan, Energy, LoadBalance, Security -> Fitness
        3 1 3 1 1 1 1;  % High Makespan, Low Energy, High LoadBalance, Low Security -> Low Fitness
        1 1 3 3 3 1 1;  % Low Makespan, Low Energy, High LoadBalance, High Security -> High Fitness
        2 2 2 2 2 1 1;  % Medium Makespan, Medium Energy, Medium LoadBalance, Medium Security -> Medium Fitness
        3 3 1 1 1 1 1;  % High Makespan, High Energy, Low LoadBalance, Low Security -> Low Fitness
        1 3 2 3 2 1 1;  % Low Makespan, High Energy, Medium LoadBalance, High Security -> Medium Fitness
        2 2 1 3 2 1 1;  % Medium Makespan, Medium Energy, Low LoadBalance, High Security -> Medium Fitness
        1 1 1 3 3 1 1;  % Low Makespan, Low Energy, Low LoadBalance, High Security -> High Fitness
        3 1 2 1 1 1 1;  % High Makespan, Low Energy, Medium LoadBalance, Low Security -> Low Fitness
        3 3 3 3 1 1 1;  % High Makespan, High Energy, High LoadBalance, High Security -> Low Fitness
        1 3 1 1 1 1 1;  % Low Makespan, High Energy, Low LoadBalance, Low Security -> Low Fitness
        2 1 3 3 3 1 1;  % Medium Makespan, Low Energy, High LoadBalance, High Security -> High Fitness
        1 2 2 1 2 1 1;  % Low Makespan, Medium Energy, Medium LoadBalance, Low Security -> Medium Fitness
        2 1 1 2 2 1 1;  % Medium Makespan, Low Energy, Low LoadBalance, Medium Security -> Medium Fitness
        1 1 2 3 3 1 1;  % Low Makespan, Low Energy, Medium LoadBalance, High Security -> High Fitness
        3 2 2 1 1 1 1;  % High Makespan, Medium Energy, Medium LoadBalance, Low Security -> Low Fitness
        2 2 3 2 3 1 1;  % Medium Makespan, Medium Energy, High LoadBalance, Medium Security -> High Fitness
        % u can Add more rules based on your complete fuzzy table
    ];
    fis = addRule(fis, rules);
    
    % Evaluate Fitness
    inputs = [makespan, energy, loadBalance, security];
    fitness = evalfis(fis, inputs);
end

% Update Positions
function positionMatrix = UpdatePositions(positionMatrix, mappedMatrix, numTasks, numVMs, fitnessValues, bestSolution)
    
    for i = 1:size(positionMatrix, 3)
        positionMatrix(:,:,i) = round(rand(numTasks, numVMs)); % Example update
    end
end

% Task Migration
function [newPosMatrix, newMappedMatrix] = MigrateTasks(positionMatrix, mappedMatrix, numTasks, numVMs)
    overloadedThreshold = 0.8; 
    underloadedThreshold = 0.2;
    overloadedVMs = sum(positionMatrix, 1) > (overloadedThreshold * numTasks / numVMs);
    underloadedVMs = sum(positionMatrix, 1) < (underloadedThreshold * numTasks / numVMs);
    for vm = 1:numVMs
        if overloadedVMs(vm)
            for task = 1:numTasks
                if positionMatrix(task, vm) == 1
                    targetVM = find(underloadedVMs, 1);
                    if ~isempty(targetVM)
                        positionMatrix(task, vm) = 0;
                        positionMatrix(task, targetVM) = 1;
                    end
                end
            end
        end
    end
    newPosMatrix = positionMatrix;
    newMappedMatrix = positionMatrix > 0;
end
function [maxObjectives, minObjectives] = UpdateMaxMinObjectives(makespan, energy, security, loadBalance, maxObjectives, minObjectives)
    maxObjectives.makespan = max(maxObjectives.makespan, makespan);
    minObjectives.makespan = min(minObjectives.makespan, makespan);
    
    maxObjectives.energy = max(maxObjectives.energy, energy);
    minObjectives.energy = min(minObjectives.energy, energy);
    
    maxObjectives.security = max(maxObjectives.security, security);
    minObjectives.security = min(minObjectives.security, security);
    
    maxObjectives.loadBalance = max(maxObjectives.loadBalance, loadBalance);
    minObjectives.loadBalance = min(minObjectives.loadBalance, loadBalance);
end

