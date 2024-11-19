
function main_scheduling()
    % Fixed parameters
    numIterations = 75;       % Number of iterations
    populationSize = 20;       % Population size
    numVMs = 20;               % Number of Virtual Machines
    taskLength = 1000;         % Task length
    
    numTasksRange = [100, 200, 500, 1000]; % Varying number of tasks
   % numTasksRange = [100]; % Varying number of tasks

    results = struct();
    
    for nTasks = numTasksRange
        fprintf('Running simulations for %d tasks...\n', nTasks);
        
        [tasks, VMs] = InitializeCharacteristics(nTasks, numVMs);
        
        % FQASO
        fprintf('  Running FQASO...\n');
        [fqasoBestSolution, fqasoBestMakespan, fqasoBestEnergy, fqasoBestSecurity, fqasoBestLoadBalance] = ...
            fqaso(tasks, VMs, populationSize, numIterations);
        results.FQASO(nTasks) = struct('BestSolution', fqasoBestSolution, ...
                                       'BestMakespan', fqasoBestMakespan, ...
                                       'BestEnergy', fqasoBestEnergy, ...
                                       'BestSecurity', fqasoBestSecurity, ...
                                       'BestLoadBalance', fqasoBestLoadBalance);
        % 
        % GAGWO
        fprintf('  Running GAGWO...\n');
        [makespanga, loadBalancega, avgSecurityLevelga, totalEnergyCostga] = gagwo(tasks, VMs, numIterations);
        results.GAGWO(nTasks) = struct('Makespan', makespanga, ...
                                       'LoadBalance', loadBalancega, ...
                                       'AvgSecurityLevel', avgSecurityLevelga, ...
                                       'TotalEnergyCost', totalEnergyCostga);
    
        %SAEA
        fprintf('  Running SAEA...\n');
        [makespansa, loadBalancesa, avgSecurityLevelsa, totalEnergyCostsa] = SAEA(tasks, VMs, numIterations);
        results.SAEA(nTasks) = struct('Makespan', makespansa, ...
                                      'LoadBalance', loadBalancesa, ...
                                      'AvgSecurityLevel', avgSecurityLevelsa, ...
                                      'TotalEnergyCost', totalEnergyCostsa);
    
        % SAEDF
        fprintf('  Running SAEDF...\n');
        [makespansae, loadBalancesae, avgSecurityLevelsae, totalEnergyCostsae] = saedf(tasks, VMs, numIterations,taskLength);
        results.SAEDF(nTasks) = struct('Makespan', makespansae, ...
                                       'LoadBalance', loadBalancesae, ...
                                       'AvgSecurityLevel', avgSecurityLevelsae, ...
                                       'TotalEnergyCost', totalEnergyCostsae);
    end
    
    fprintf('\nSummary of Results:\n');
    disp(results);
end

% Initialize Tasks and VM Characteristics
function [tasks, VMs] = InitializeCharacteristics(numTasks, numVMs)
    tasks = struct();
    tasks.sensitivity.confidentiality = rand(1, numTasks);  
    tasks.sensitivity.integrity = rand(1, numTasks);        
    tasks.sensitivity.authentication = rand(1, numTasks);   
    tasks.executionTime = randi([50, 200], 1, numTasks);    
    tasks.resource.CPU = randi([1, 10], 1, numTasks);       
    tasks.resource.memory = randi([1, 16], 1, numTasks);    
    
    % VM Characteristics
    VMs = struct();
    VMs.trustLevel = rand(1, numVMs);                      
    VMs.capacity.CPU = randi([10, 100], 1, numVMs);       
    VMs.capacity.memory = randi([16, 64], 1, numVMs);     
    VMs.power.idle = randi([40, 100], 1, numVMs);          
    VMs.power.busy = randi([150, 300], 1, numVMs);         
end