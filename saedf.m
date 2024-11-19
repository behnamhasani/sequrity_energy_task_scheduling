
function [makespan, loadBalance, avgSecurityLevel, totalEnergyCost] = saedf(tasks, VMs, timeSlots,taskLength)
    numTasks = length(tasks.executionTime);
    numServers = length(VMs.capacity.CPU);
    tasks.dataSizes = taskLength * ones(1, numTasks); % Ensure dataSizes is fully initialized
    queueBacklog = zeros(1, timeSlots + 1); % Initialize backlog queue
    queueBacklog(1) = 50;                  % Initial workload backlog
    totalEnergyCost = 0;                   %  energy cost
    totalSecurityMetric = 0;               %  security achieved
    serverLoads = zeros(1, numServers);    % Load per server
    electricityPrices = rand(1, timeSlots) * 0.5 + 0.2; %  electricity prices
    V = 10;                                % Lyapunov control parameter
    f_min = 1.2; f_max = 3.2;              % Frequency bounds
    alpha = 6.1; delta = 100; theta = 3;   % Energy parameters

    % Main simulation loop
    for t = 1:timeSlots
        fprintf('Time Slot %d/%d\n', t, timeSlots);
        
        p_t = electricityPrices(t);
        
        %  new workload (Poisson process)
        jobsArriving = poissrnd(10);
        jobWorkloads = randi([10, 50], 1, jobsArriving);
        
        %  security workload for new jobs
        [securityWorkloads, avgSecurityLevels] = CalculateSecurityWorkload(tasks, jobsArriving);
        
        %  including security overhead
        totalWorkload = sum(jobWorkloads .* (1 + securityWorkloads(1:jobsArriving)));
        
        [f_opt, energyCost] = OptimizeEnergy(totalWorkload, numServers, p_t, V, f_min, f_max, alpha, delta, theta);
        totalEnergyCost = totalEnergyCost + energyCost;
        
        % Distribute workload across servers
        serverLoads = DistributeWorkload(serverLoads, totalWorkload, numServers);
        
        % Update backlog queue
        serviceRate = f_opt * numServers;
        queueBacklog(t + 1) = max(queueBacklog(t) - serviceRate, 0) + totalWorkload;
        
        % Accumulate security metrics
        totalSecurityMetric = totalSecurityMetric + mean(avgSecurityLevels);
    end

    makespan = max(serverLoads);
    avgLoad = mean(serverLoads);
    loadBalance = 1 - (sqrt(sum((serverLoads - avgLoad).^2) / numServers) / avgLoad);
    avgSecurityLevel = totalSecurityMetric / timeSlots;

    fprintf('Total Energy Cost: %.2f\n', totalEnergyCost);
    fprintf('Final Queue Backlog: %.2f\n', queueBacklog(end));
    fprintf('Makespan: %.2f\n', makespan);
    fprintf('Load Balance: %.2f\n', loadBalance);
    fprintf('Average Security Level: %.2f\n', avgSecurityLevel);
end


function [securityWorkloads, avgSecurityLevels] = CalculateSecurityWorkload(tasks, jobsArriving)
    % Security parameters
    securityLevelsAuth = [0, 0.55, 0.91, 1.0];
    securityLevelsInt = [0, 0.18, 0.26, 0.36, 0.45, 0.63, 1.0];
    securityLevelsConf = [0, 0.08, 0.14, 0.36, 0.46, 0.64, 1.0];
    betaAuth = 1600; betaInt = 2400; betaConf = 800;

    securityWorkloads = zeros(1, jobsArriving);
    avgSecurityLevels = zeros(1, jobsArriving);

    for j = 1:jobsArriving
        sla = securityLevelsAuth(randi([1, length(securityLevelsAuth)]));
        slg = securityLevelsInt(randi([1, length(securityLevelsInt)]));
        slc = securityLevelsConf(randi([1, length(securityLevelsConf)]));

        Din = tasks.dataSizes(j);
        Dout = tasks.dataSizes(j);
        SW_auth = betaAuth * sla;
        SW_int = betaInt * slg * Din;
        SW_conf = betaConf * slc * Dout;

        securityWorkloads(j) = SW_auth + SW_int + SW_conf;
        avgSecurityLevels(j) = mean([sla, slg, slc]);
    end
end

% Optimize Energy
function [f_opt, energyCost] = OptimizeEnergy(totalWorkload, numServers, p_t, V, f_min, f_max, alpha, delta, theta)
    f_opt = max(min(V * sqrt(totalWorkload / numServers), f_max), f_min);
    powerConsumption = alpha * f_opt^theta + delta;
    energyCost = powerConsumption * p_t;
end

% Distribute Workload
function serverLoads = DistributeWorkload(serverLoads, totalWorkload, numServers)
    workloadPerServer = totalWorkload / numServers;
    serverLoads = serverLoads + workloadPerServer;
end
