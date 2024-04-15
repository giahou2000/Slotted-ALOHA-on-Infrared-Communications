qLow = 0.1; % probability
PLow = 0.0001; % (w)
RLow = 0.001; % (bps/Hz)
qHigh = 1; % probability
PHigh = 0.05; % (W)
RHigh = 2.5; % (bps/Hz)
blo = [qLow, PLow, RLow]; % minimum limit for the tested variables
bup = [qHigh, PHigh, RHigh]; % maximum limit for the tested variables
velocity_low = -abs(bup - blo);
velocity_high = abs(bup - blo);
num_devices = 5; % Change this to the desired number of devices (K)
num_variables = 3;
num_particles = 50; % number of particles that will search for the best position
phi_p = 1.1; % cognitive parameter
phi_g = 1.1; % social parameter
w = 0.7; % inertia weight
global_best_fitness = 0;
s = 1e-08;
heta = 0.6;
H_0 = [1.6749e-06, 2.0559e-06, 4.6345e-06, 2.4155e-06, 1.6482e-06];

% Initialize particles
% Pre-allocate empty arrays of structs
particles = struct.empty(num_particles, 0);  % 0 indicates no pre-defined fields
global_best = struct.empty(num_devices, 0);

p = 0;
for max_iterations = 100:100:5000
    p = p+1;
end
converge = struct.empty(p, 0);


% Iterations for convergance
p = 1;
for max_iterations = 100:100:5000
    
    % Access and modify individual structs
    for i = 1:num_particles
        for j = 1:num_devices
          particles(i).node(j).position(1) = qLow + (qHigh - qLow) * rand; % Random initial probability
          particles(i).node(j).position(2) = PLow + (PHigh - PLow) * rand; % Random initial power
          particles(i).node(j).position(3) = RLow + (RHigh - RLow) * rand; % Random initial transmission rate
          particles(i).node(j).velocity = velocity_low + (velocity_high - velocity_low) * rand; % Random initial velocities
          particles(i).node(j).best_position = particles(i).node(j).position; % Best known positions
        end
        particles(i).best_fitness = fitness(particles, num_devices, i, s, heta, H_0);
    end
    
    % Initialize global best position and fitness
    for i = 1:num_particles
        if particles(i).best_fitness > global_best_fitness
            global_best_fitness = particles(i).best_fitness;
            for j = 1:num_devices
                global_best(j).position = particles(i).node(j).position;
            end
        end
    end
    
    
    
    for iteration = 1:max_iterations
    
        % For each particle
        for i = 1:num_particles
    
            % For each node of the particle
            for j = 1:num_devices
                % Update velocity
                rp = rand(1, num_variables);
                rg = rand(1, num_variables);
                old_velocity = particles(i).node(j).velocity;
                particles(i).node(j).velocity = w * old_velocity + phi_p * rp .* (particles(i).node(j).best_position - particles(i).node(j).position) + phi_g * rg .* (global_best(j).position - particles(i).node(j).position);
                old_position = particles(i).node(j).position;
                temp = old_position + particles(i).node(j).velocity;
                if (temp(1) > qLow) && (temp(1) < qHigh) && (temp(2) > PLow) && (temp(2) < PHigh) && (temp(3) > RLow) && (temp(3) < RHigh)% implement the limits of the variables in if statement
                    particles(i).node(j).position = temp;
                else
                    continue
                end
            end
            
            % Evaluate fitnessprint
            current_fitness = fitness(particles, num_devices, i, s, heta, H_0);
            
            % Update personal best
            if current_fitness > particles(i).best_fitness
                particles(i).best_fitness = current_fitness;
                for j = 1:num_devices
                    particles(i).node(j).best_position = particles(i).node(j).position;
                end
                if current_fitness > global_best_fitness
                    global_best_fitness = particles(i).best_fitness;
                    for j = 1:num_devices
                        global_best(j).position = particles(i).node(j).position;
                    end
                end
            end
        end
        
        % Display current best fitness value for each iteration
        fprintf('Iteration %d: Current best fitness = %.4f\n', iteration, global_best_fitness);
    end
    

    converge(p).global_best = global_best;
    converge(p).best_fitness = global_best_fitness;
    disp("Saved results!!!!!!!!!!!!")
    p = p + 1;
end








%% Supporting Equations for calculating the objective function
% ______________________________________________

% Channel Statistical Model (Gamma distribution)
% To describe the channel DC gain we use the gamma distribution
% x = 1:10;
% a = 13.79;
% b = 0.04;
% f_x = gampdf(x, a, b);

% Objective function
function value = fitness(particles, nodes_num, which_particle, s, heta, H)
    % value = rand; % just for testing the functionality of the rest of the
    % algorithm
    % Compute the Rk_hut and the Pk using the functions below
    value = 0;
    for i = 1:nodes_num
        Rk_power = avRate (particles(which_particle).node(i).position(3), i, particles(which_particle).node(i).position(2), s, heta, H(i));
        Rk_hat = avThrouput (i, Rk_power, particles(which_particle));
        Pk = particles(which_particle).node(i).position(2);
        value = value + (Rk_hat/Pk);
    end
end

% Average throughput of the network Rk_hat
% k is the number of the kth node
% Rk_power is the average rate of the node
% q(k) is the probability of channel access of the kth node
function Rk_hat = avThrouput (k, Rk_power, particle)
    temp1 = Rk_power * particle.node(k).position(1);
    temp2 = 1;
    for i = 1:k
        if i ~= k
            temp2 = temp2 * (1 - particle.node(i).position(1));
        end
    end
    Rk_hat = temp1 * temp2;
end

% Average throughput of the network Rk_hat but for the results calculations
% We do that so that we use the global_best struct instead of the particle
% k is the number of the kth node
% Rk_power is the average rate of the node
% q(k) is the probability of channel access of the kth node
function Rk_hat = avThrouputforSNR (k, Rk_power, particle)
    temp1 = Rk_power * particle(k).position(1);
    temp2 = 1;
    for i = 1:k
        if i ~= k
            temp2 = temp2 * (1 - particle(i).position(1));
        end
    end
    Rk_hat = temp1 * temp2;
end


% Average rate of the kth node Rk_power

function Rk_power = avRate (Rk, ~, Pk, s, heta, h0)
    Xk = Xk_helper(s, Rk, Pk, heta, h0);
    a = 6.11;
    b = 0.07;
    g = gammainc(Xk/(b^2), a);
    % G = gamma(k);
    % temp1 = g/G;
    % temp2 = 1 - temp1;
    temp2 = 1 - g;
    Rk_power = Rk * temp2;
end


% Xk helper function

function y = Xk_helper(s, Rk, Pk, heta, h0)
    numerator = 2 * pi * (s^2) * ((2 ^ Rk) - 1);
    denominator = exp(1) * ((abs(heta * h0 * Pk))^2);
    fraction = numerator/denominator;
    if fraction < 0
        y = 0;
    else
        y = sqrt(fraction);
    end
end


