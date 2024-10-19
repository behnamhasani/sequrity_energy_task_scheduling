

numServers = 10;        % Number of homogeneous servers in IDC
timeSlots = 100;        % Total time slots for simulation
electricityPrices = rand(1, timeSlots); % Random electricity prices per time slot
queueBacklog = zeros(1, timeSlots);     % Queue backlog to store workload to be processed

% Security Levels and Parameters 
securityLevelsAuth = [0, 0.55, 0.91, 1.0];        % Authentication levels
securityLevelsInt = [0, 0.18, 0.26, 0.36, 0.45, 0.63, 1.0]; % Integrity levels
securityLevelsConf = [0, 0.08, 0.14, 0.36, 0.46, 0.64, 1.0]; % Confidentiality levels

% Constants for security workload calculations
betaAuth = 1600;  % Authentication constant
betaInt = 2400;   % Integrity constant
betaConf = 800;   % Confidentiality constant

% Task Parameters
numTasks = 100;       % Number of tasks
dataSizes = rand(1, numTasks);  % Data sizes (in MI) for tasks

% Control parameters for SAEDF (Lyapunov optimization)
V = 10; % Lyapunov control parameter
Pmax = 1000; % Max power consumption
f_min = 1.2; % Min server frequency
f_max = 3.2; % Max server frequency
theta = 3;   % Power function exponent
alpha = 6.1; % Proportionality constant for energy
delta = 100; % Idle power consumption



% Initialize queue and workload arrays
queueBacklog(1) = 50;  % Initial workload
totalEnergyCost = 0;   % Energy cost for IDC
totalSecurityWorkload = zeros(1, numTasks);  % Security workload for each task

% Initialize random security levels for tasks
taskAuthLevels = randi([1, length(securityLevelsAuth)], 1, numTasks); % Authentication levels
taskIntLevels = randi([1, length(securityLevelsInt)], 1, numTasks);   % Integrity levels
taskConfLevels = randi([1, length(securityLevelsConf)], 1, numTasks); % Confidentiality levels

for t = 1:timeSlots
    % Step 1: Get the electricity price for the current time slot
    p_t = electricityPrices(t);
    
    % Step 2: Workload Arrival - Poisson process (Random number of jobs arrive)
    jobsArriving = poissrnd(10); % Poisson distributed job arrivals
    jobWorkloads = randi([10, 50], 1, jobsArriving); % Random job workloads
    jobSecurityLevels = randi([1, 3], 1, jobsArriving); % Random security levels for jobs

 
    for j = 1:jobsArriving
        % Get the security levels for authentication, integrity, and confidentiality
        sla = securityLevelsAuth(taskAuthLevels(j));  % Authentication level
        slg = securityLevelsInt(taskIntLevels(j));    % Integrity level
        slc = securityLevelsConf(taskConfLevels(j));  % Confidentiality level
        
        % Get the size of data to be protected
        Din = dataSizes(j);  % Input data size
        Dout = dataSizes(j); % Output data size

        % Calculate security workload for each service
        SW_auth = betaAuth * sla;  % Authentication workload
        SW_int = betaInt * slg * Din;  % Integrity workload (based on input data)
        SW_conf = betaConf * slc * Dout; % Confidentiality workload (based on output data)

        % Total security workload for the task
        totalSecurityWorkload(j) = SW_auth + SW_int + SW_conf;
    end

    % Step 3: Total Workload including security overhead
    totalWorkload = sum(jobWorkloads .* (1 + totalSecurityWorkload(1:jobsArriving)));
    
    % Update queue backlog
    queueBacklog(t+1) = max(queueBacklog(t) - totalWorkload, 0) + totalWorkload;

    % Calculate required resource based on the total workload including security overhead
    resourceRequired = totalWorkload; % Simplified to total workload including security


    % Apply Dynamic Voltage and Frequency Scaling (DVFS) to adjust server frequency
    % The objective is to choose an optimal frequency that minimizes the cost function.

    f_opt = max(min(V * sqrt(resourceRequired), f_max), f_min); % Optimized frequency based on Lyapunov drift
    powerConsumption = alpha * f_opt^theta + delta; % Power function model
    energyCost = powerConsumption * p_t; % Energy cost at this time slot based on the power consumption and electricity price

    % Update queue backlog
    % Calculate the effective service rate, i.e., how much workload is processed
    serviceRate = f_opt * numServers; % The total processing capacity of the IDC at optimal frequency
    queueBacklog(t+1) = max(queueBacklog(t) - serviceRate, 0) + totalWorkload; % Update the queue backlog

    % Accumulate the total energy cost over all time slots
    totalEnergyCost = totalEnergyCost + energyCost;
end

% ============================
% Output the Results
% ============================
disp(['Total Energy Cost: ', num2str(totalEnergyCost)]);
disp(['Final Queue Backlog: ', num2str(queueBacklog(end))]);
disp('Total Security Workload for each task:');
disp(totalSecurityWorkload);

