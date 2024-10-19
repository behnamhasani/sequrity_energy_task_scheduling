function [bestSolution, bestMakespan, bestEnergy, bestSecurity, bestLoadBalance, iterationResults] = fqaso(numTasks, numVMs, populationSize, maxIterations)
    % Initialize task and VM characteristics
    [tasks, VMs] = InitializeCharacteristics(numTasks, numVMs);
    
    % Initialize population and other parameters
    positionMatrix = InitializePopulation(numTasks, numVMs, populationSize);
    mappedMatrix = MapMatrix(positionMatrix, numTasks, numVMs);
    
    % Initialize best solution and objectives placeholders
    bestSolution = [];
    bestFitness = Inf;
    
    bestMakespan = Inf;
    bestEnergy = Inf;
    bestSecurity = Inf;
    bestLoadBalance = Inf;

        % Initialize maximum and minimum values for the objectives
    maxObjectives.makespan = -inf;
    maxObjectives.energy = -inf;
    maxObjectives.security = -inf;
    maxObjectives.loadBalance = -inf;
    
    minObjectives.makespan = inf;
    minObjectives.energy = inf;
    minObjectives.security = inf;
    minObjectives.loadBalance = inf;
    
    % Array to store results at each iteration
    iterationResults = struct('iteration', [], 'bestSolution', [], 'bestMakespan', [], 'bestEnergy', [], 'bestSecurity', [], 'bestLoadBalance', []);
    
    % Begin iterations
    for t = 1:maxIterations
        for i = 1:populationSize
            % Evaluate fitness for each q-atom
            [fitness, makespan, energy, security, loadBalance] = EvaluateFitness(positionMatrix(:,:,i), tasks, VMs,maxObjectives,minObjectives);
            fitnessValues(i) = fitness;

            % Update the best solution if the current one is better
            if fitness < bestFitness
                bestFitness = fitness;
                bestSolution = positionMatrix(:,:,i);
                bestMakespan = makespan;
                bestEnergy = energy;
                bestSecurity = security;
                bestLoadBalance = loadBalance;
            end
        end
        
        % Store the best values of this iteration
        iterationResults(t).iteration = t;
        iterationResults(t).bestSolution = bestSolution;
        iterationResults(t).bestMakespan = bestMakespan;
        iterationResults(t).bestEnergy = bestEnergy;
        iterationResults(t).bestSecurity = bestSecurity;
        iterationResults(t).bestLoadBalance = bestLoadBalance;
        
        % Update positions using ASO mechanism
        positionMatrix = UpdatePositions(positionMatrix, mappedMatrix, numTasks, numVMs, fitnessValues, bestSolution);
        
        % Perform task migration for load balancing
        [positionMatrix, mappedMatrix] = MigrateTasks(positionMatrix, mappedMatrix, numTasks, numVMs);
    end
    
    % Output the best results
    disp('Best solution found:');
    disp(bestSolution);
    disp(['Makespan: ', num2str(bestMakespan)]);
    disp(['Energy Consumption: ', num2str(bestEnergy)]);
    disp(['Security Bias: ', num2str(bestSecurity)]);
    disp(['Load Balance: ', num2str(bestLoadBalance)]);
end


function [tasks, VMs] = InitializeCharacteristics(numTasks, numVMs)
    % Initialize tasks' characteristics
    tasks = struct();
    tasks.sensitivity.confidentiality = rand(1, numTasks);  % Confidentiality (0 to 1)
    tasks.sensitivity.integrity = rand(1, numTasks);        % Integrity (0 to 1)
    tasks.sensitivity.authentication = rand(1, numTasks);   % Authentication (0 to 1)
    tasks.executionTime = randi([50, 200], 1, numTasks);    % Task execution time (random range)
    tasks.resource.CPU = randi([1, 10], 1, numTasks);       % CPU requirement (random range)
    tasks.resource.memory = randi([1, 16], 1, numTasks);    % Memory requirement in GB (random)

    % Initialize VMs' characteristics
    VMs = struct();
    VMs.trustLevel = rand(1, numVMs);                      % Trust level (0 to 1)
    VMs.capacity.CPU = randi([10, 100], 1, numVMs);        % CPU capacity of each VM
    VMs.capacity.memory = randi([16, 64], 1, numVMs);      % Memory capacity in GB
    VMs.power.idle = randi([40, 100], 1, numVMs);          % Idle power consumption (in watts)
    VMs.power.busy = randi([150, 300], 1, numVMs);         % Busy power consumption (in watts)
end

% Initialize the population of q-atoms (solutions)
function positionMatrix = InitializePopulation(numTasks, numVMs, populationSize)
    % Initialize random positions for each task-VM combination for the entire population
    positionMatrix = randi([0 1], numTasks, numVMs, populationSize); % Random binary allocation of tasks to VMs
end


% Map matrix for q-atoms (binary map indicating task allocation to VMs)
function mappedMatrix = MapMatrix(positionMatrix, numTasks, numVMs)
    mappedMatrix = positionMatrix > 0; % Binary map
end

% Evaluate the fitness of a solution
function [fitness, makespan, energy, security, loadBalance] = EvaluateFitness(positionMatrix, tasks, VMs,maxObjectives, minObjectives)
    
    makespan = CalculateMakespan(positionMatrix, tasks, VMs);
    
  
    security = CalculateSecurityBias(positionMatrix, tasks, VMs);
    
    
    energy = CalculateEnergyConsumption(positionMatrix, tasks, VMs);
    
    
    loadBalance = CalculateLoadBalance(positionMatrix, tasks, VMs);
    

    [maxObjectives, minObjectives] = UpdateMaxMinObjectives(makespan, energy, security, loadBalance, maxObjectives, minObjectives);

    % Combine objectives into a single fitness value using fuzzy system
    fitness = FuzzySystem(makespan, energy, loadBalance, security,maxObjectives,minObjectives);
end

function fitness = FuzzySystem(makespan, energy, loadBalance, security,maxObjectives,minObjectives)
    % Create fuzzy inference system (FIS) object
    fis = mamfis('Name', 'FitnessFIS');
    
    % Add input variables
    fis = addInput(fis, [0 1], 'Name', 'Makespan');  % Normalize all inputs between 0 and 1
    fis = addInput(fis, [0 1], 'Name', 'Energy');
    fis = addInput(fis, [0 1], 'Name', 'LoadBalance');
    fis = addInput(fis, [0 1], 'Name', 'Security');
    
    % Add output variable for fitness
    fis = addOutput(fis, [0 1], 'Name', 'Fitness');
    
    % Define fuzzy membership functions for input variables
    fis = addMF(fis, 'Makespan', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Makespan', 'trapmf', [0.4 0.6 1 1], 'Name', 'High');
    
    fis = addMF(fis, 'Energy', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Energy', 'trapmf', [0.4 0.6 1 1], 'Name', 'High');
    
    fis = addMF(fis, 'LoadBalance', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'LoadBalance', 'trapmf', [0.4 0.6 1 1], 'Name', 'High');
    
    fis = addMF(fis, 'Security', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Security', 'trapmf', [0.4 0.6 1 1], 'Name', 'High');
    
    % Define fuzzy membership functions for output variable
    fis = addMF(fis, 'Fitness', 'trapmf', [0 0 0.3 0.5], 'Name', 'Low');
    fis = addMF(fis, 'Fitness', 'trapmf', [0.4 0.6 1 1], 'Name', 'High');
    
    % Define fuzzy rules (correct length for m+n+2)
    ruleList = [
        1 1 1 1 1 1 1;  % If makespan is Low, energy is Low, load balance is Low, and security is Low -> Fitness is Low
        1 1 1 2 1 1 1;  % If makespan is Low, energy is Low, load balance is Low, and security is High -> Fitness is Low
        1 1 2 1 1 1 1;  % If makespan is Low, energy is Low, load balance is High, and security is Low -> Fitness is Low
        1 1 2 2 2 1 1;  % If makespan is Low, energy is Low, load balance is High, and security is High -> Fitness is Medium
        1 2 1 1 1 1 1;  % If makespan is Low, energy is High, load balance is Low, and security is Low -> Fitness is Low
        1 2 1 2 2 1 1;  % If makespan is Low, energy is High, load balance is Low, and security is High -> Fitness is Medium
        1 2 2 1 2 1 1;  % If makespan is Low, energy is High, load balance is High, and security is Low -> Fitness is Medium
        1 2 2 2 2 1 1;  % If makespan is Low, energy is High, load balance is High, and security is High -> Fitness is High
        2 1 1 1 1 1 1;  % If makespan is High, energy is Low, load balance is Low, and security is Low -> Fitness is Low
        2 1 1 2 2 1 1;  % If makespan is High, energy is Low, load balance is Low, and security is High -> Fitness is Medium
        2 1 2 1 2 1 1;  % If makespan is High, energy is Low, load balance is High, and security is Low -> Fitness is Medium
        2 1 2 2 2 1 1;  % If makespan is High, energy is Low, load balance is High, and security is High -> Fitness is Medium
        2 2 1 1 2 1 1;  % If makespan is High, energy is High, load balance is Low, and security is Low -> Fitness is Medium
        2 2 1 2 2 1 1;  % If makespan is High, energy is High, load balance is Low, and security is High -> Fitness is Medium
        2 2 2 1 2 1 1;  % If makespan is High, energy is High, load balance is High, and security is Low -> Fitness is Medium
        2 2 2 2 2 1 1;  % If makespan is High, energy is High, load balance is High, and security is High -> Fitness is High
    ];
    
    % Apply rules to FIS
    fis = addRule(fis, ruleList);
    
 normalizedMakespan = NormalizeObjective(makespan, minObjectives.makespan, maxObjectives.makespan);
    normalizedEnergy = NormalizeObjective(energy, minObjectives.energy, maxObjectives.energy);
    normalizedLoadBalance = NormalizeObjective(loadBalance, minObjectives.loadBalance, maxObjectives.loadBalance);
    normalizedSecurity = NormalizeObjective(security, minObjectives.security, maxObjectives.security);

    % Evaluate fitness based on the normalized inputs
    inputs = [normalizedMakespan, normalizedEnergy, normalizedLoadBalance, normalizedSecurity];
    fitness = evalfis(fis, inputs);
end

% Function to normalize inputs to [0, 1] range
function normalizedValue = NormalizeInput(value, minVal, maxVal)
    
    normalizedValue = (value - minVal) / (maxVal - minVal);
    normalizedValue = min(max(normalizedValue, 0), 1);  % Ensure the value is within [0, 1]
end




function makespan = CalculateMakespan(positionMatrix, tasks, VMs)
    vmCompletionTimes = zeros(1, length(VMs.trustLevel));  
    
    % Loop through each VM to calculate task completion times
    for vm = 1:length(VMs.trustLevel)
        assignedTasks = find(positionMatrix(:, vm) == 1); % Tasks assigned to this VM
        
        if ~isempty(assignedTasks)
            
            for task = assignedTasks'
                
                % A more powerful VM (higher CPU capacity) reduces task execution time
                adjustedExecutionTime = tasks.executionTime(task) * (tasks.resource.CPU(task) / VMs.capacity.CPU(vm));
                
                % Accumulate the total completion time for the VM
                vmCompletionTimes(vm) = vmCompletionTimes(vm) + adjustedExecutionTime;
            end
        end
    end
    
    
    makespan = max(vmCompletionTimes); % The makespan is the maximum VM completion time
end

% Calculate energy consumption based on task execution and VM power usage
function energy = CalculateEnergyConsumption(positionMatrix, tasks, VMs)
    energy = 0;
    
    % Loop through each VM
    for vm = 1:length(VMs.trustLevel)
        assignedTasks = find(positionMatrix(:, vm) == 1); % Tasks assigned to this VM
        
        if ~isempty(assignedTasks)
            % If tasks are assigned, VM is busy
            T_busy = sum(tasks.executionTime(assignedTasks)); % Busy time is sum of task execution times
            energy_busy = VMs.power.busy(vm) * T_busy;       % Busy energy
            
            % Security energy overhead
            E_security = sum(tasks.sensitivity.confidentiality(assignedTasks)) + ...
                         sum(tasks.sensitivity.integrity(assignedTasks)) + ...
                         sum(tasks.sensitivity.authentication(assignedTasks));
                     
            % Add busy energy and security overhead
            energy = energy + energy_busy + E_security;
        else
            % VM is idle
            T_idle = 200; % Example idle time
            energy_idle = VMs.power.idle(vm) * T_idle; % Idle energy
            energy = energy + energy_idle;
        end
    end
end

% Calculate security bias based on task sensitivity and VM trust levels
function security = CalculateSecurityBias(positionMatrix, tasks, VMs)
    security = 0;
    
    % Loop through each task and calculate security bias
    for task = 1:length(tasks.executionTime)
        for vm = 1:length(VMs.trustLevel)
            if positionMatrix(task, vm) == 1
                confidentialityOverhead = tasks.sensitivity.confidentiality(task) * (1 - VMs.trustLevel(vm));
                integrityOverhead = tasks.sensitivity.integrity(task) * (1 - VMs.trustLevel(vm));
                authenticationOverhead = tasks.sensitivity.authentication(task) * (1 - VMs.trustLevel(vm));
                
                % Accumulate total security bias
                security = security + (confidentialityOverhead + integrityOverhead + authenticationOverhead);
            end
        end
    end
end

% Calculate load balance based on task distribution across VMs
function loadBalance = CalculateLoadBalance(positionMatrix, tasks, VMs)
    % Calculate the resource load (CPU and memory) on each VM
    vmLoad.CPU = zeros(1, length(VMs.trustLevel));
    vmLoad.memory = zeros(1, length(VMs.trustLevel));
    
    % Loop through each VM
    for vm = 1:length(VMs.trustLevel)
        assignedTasks = find(positionMatrix(:, vm) == 1); % Find tasks assigned to this VM
        
        % Calculate the total CPU and memory load on the VM from assigned tasks
        vmLoad.CPU(vm) = sum(tasks.resource.CPU(assignedTasks));        % Total CPU load from assigned tasks
        vmLoad.memory(vm) = sum(tasks.resource.memory(assignedTasks));  % Total memory load from assigned tasks
    end
    
    
    cpuUtilization = vmLoad.CPU ./ VMs.capacity.CPU;         % CPU utilization as fraction of VM capacity
    memoryUtilization = vmLoad.memory ./ VMs.capacity.memory; % Memory utilization as fraction of VM capacity
    
    % Average utilization across all VMs
    avgCpuUtilization = mean(cpuUtilization);
    avgMemoryUtilization = mean(memoryUtilization);
    
    % Calculate load imbalance as the squared deviation from average utilization
    cpuImbalance = sum((cpuUtilization - avgCpuUtilization).^2);
    memoryImbalance = sum((memoryUtilization - avgMemoryUtilization).^2);
    
    % Combine CPU and memory imbalances to get total load imbalance
    totalImbalance = cpuImbalance + memoryImbalance;
    
    loadBalance = 1 / (1 + totalImbalance);  % Adding 1 to avoid division by zero
end


function positionMatrix = UpdatePositions(positionMatrix, mappedMatrix, numTasks, numVMs, fitnessValues, bestSolution)
    
    
   
    c1 = 1.5;  % Attraction coefficient
    c2 = 2.0;  % Repulsion coefficient
    alpha = 0.2;  % A constant value for the velocity update
    beta = 0.5;  % A constant value for position update
    atomMass = fitnessValues;  % Mass of each atom is proportional to its fitness

    % Initialize variables
    numAtoms = size(positionMatrix, 3);  % Number of atoms (population size)
    velocityMatrix = zeros(size(positionMatrix));  % Velocity for each solution

    % Calculate forces between atoms (solutions)
    for i = 1:numAtoms
        forceAttraction = zeros(numTasks, numVMs);
        forceRepulsion = zeros(numTasks, numVMs);

        % Compute attraction and repulsion forces between atoms
        for j = 1:numAtoms
            if i ~= j  % Skip the same atom
                distance = norm(positionMatrix(:,:,i) - positionMatrix(:,:,j), 2);  % Euclidean distance

                % Attractive force (based on fitness - better solutions attract others)
                attractiveForce = (atomMass(j) / distance^2) * (positionMatrix(:,:,j) - positionMatrix(:,:,i));
                forceAttraction = forceAttraction + c1 * attractiveForce;

                % Repulsive force (to avoid stagnation - keeps solutions diverse)
                repulsiveForce = (1 / distance^2) * (positionMatrix(:,:,i) - positionMatrix(:,:,j));
                forceRepulsion = forceRepulsion + c2 * repulsiveForce;
            end
        end

        % Update velocity for each atom (solution)
        velocityMatrix(:,:,i) = alpha * velocityMatrix(:,:,i) + beta * (forceAttraction + forceRepulsion);
        
        % Update position based on the new velocity
        newPosition = positionMatrix(:,:,i) + velocityMatrix(:,:,i);
        
        % Ensure the positions are valid (0 or 1 for binary assignments)
        positionMatrix(:,:,i) = round(min(max(newPosition, 0), 1));
    end
end




function normalizedValue = NormalizeObjective(value, minVal, maxVal)
    if maxVal == minVal
        normalizedValue = 0.5;  % If no variation, default to 0.5
    else
        normalizedValue = (value - minVal) / (maxVal - minVal);
    end
    normalizedValue = min(max(normalizedValue, 0), 1);  % Ensure within [0, 1]
end


function rotation_angle = CalculateRotationAngle(currentPosition, bestPosition, alpha, beta, t, maxIterations)
    % Define the coordinate rotation angle (Eq. 47)
    rotation_angle = alpha * (bestPosition - currentPosition) + beta * (rand(size(currentPosition)) - 0.5);
    
    % Scale the rotation angle based on iteration (Eq. 48)
    scale_factor = (1 - t / maxIterations);
    rotation_angle = rotation_angle * scale_factor;
end

% Task migration function for load balancing
function [newPosMatrix, newMappedMatrix] = MigrateTasks(positionMatrix, mappedMatrix, numTasks, numVMs)
    % Identify overloaded and underloaded VMs for task migration
    overloadedThreshold = 0.8; % Threshold for overloaded VMs
    underloadedThreshold = 0.2; % Threshold for underloaded VMs
    
    overloadedVMs = sum(positionMatrix, 1) > (overloadedThreshold * numTasks / numVMs); % Find overloaded VMs
    underloadedVMs = sum(positionMatrix, 1) < (underloadedThreshold * numTasks / numVMs); % Find underloaded VMs
    
    % Perform simple migration from overloaded to underloaded VMs
    for vm = 1:numVMs
        if overloadedVMs(vm)
            % Migrate tasks to underloaded VMs
            for task = 1:numTasks
                if positionMatrix(task, vm) == 1
                    targetVM = find(underloadedVMs, 1); % Find the first underloaded VM
                    if ~isempty(targetVM)
                        positionMatrix(task, vm) = 0; % Remove task from overloaded VM
                        positionMatrix(task, targetVM) = 1; % Assign task to underloaded VM
                        underloadedVMs(targetVM) = false; % Mark the target VM as loaded
                    end
                end
            end
        end
    end
    
    % Return the updated position matrix and mapping matrix
    newPosMatrix = positionMatrix;
    newMappedMatrix = positionMatrix > 0; % Update the mapped matrix
end

% Function to update the maximum and minimum objectives
function [maxObjectives, minObjectives] = UpdateMaxMinObjectives(makespan, energy, security, loadBalance, maxObjectives, minObjectives)
    % Objective 1: Makespan
    maxObjectives.makespan = max(maxObjectives.makespan, makespan);
    minObjectives.makespan = min(minObjectives.makespan, makespan);
    
    % Objective 2: Energy Consumption
    maxObjectives.energy = max(maxObjectives.energy, energy);
    minObjectives.energy = min(minObjectives.energy, energy);
    
    % Objective 3: Security Bias
    maxObjectives.security = max(maxObjectives.security, security);
    minObjectives.security = min(minObjectives.security, security);
    
    % Objective 4: Load Balance
    maxObjectives.loadBalance = max(maxObjectives.loadBalance, loadBalance);
    minObjectives.loadBalance = min(minObjectives.loadBalance, loadBalance);
end
