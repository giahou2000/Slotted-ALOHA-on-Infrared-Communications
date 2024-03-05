%% Performance of Slotted ALOHA on infrared WBANs
% This simulation tries to describe the performance of Slotted ALOHA
% when applied as a MAC protocol for a WBAN system that consists of
% 6 nodes. For the communication medium it is assumed that we use the
% infrared part of the electromagnetic spectrum.

% Authors:
% Giachoudis Christos, christos.giachoudis@fresnel.fr
% Vasilis Papanikolaou,
% Konstantinos Rallis,



%% PSO - Particle Swarm Optimization
% A brute force / exploration method for the optimization problem

% PSO Parameters
num_devices = 5; % Change this to the desired number of devices (K)
num_dimensions = 3 * num_devices;
num_particles = 30;
max_iterations = 100;
c1 = 2; % cognitive parameter
c2 = 2; % social parameter
inertia_weight = 0.7;

% Initialize particles
particles.position = rand(num_particles, num_dimensions); % Random initial positions
particles.velocity = rand(num_particles, num_dimensions); % Random initial velocities
particles.best_position = particles.position; % Best known positions
% particles.best_fitness = zeros(num_particles, 1); % Best known fitness values

% Evaluate fitness for each particle
for i = 1:num_particles
    particles.best_fitness(i) = your_objective_function(particles.best_position(i, :), num_devices);
end

% Initialize global best position and fitness
[global_best_fitness, global_best_index] = min(particles.best_fitness);
global_best_position = particles.best_position(global_best_index, :);

% Main PSO loop
for iteration = 1:max_iterations
    % Update particle velocities and positions
    for i = 1:num_particles
        r1 = rand(1, num_dimensions);
        r2 = rand(1, num_dimensions);
        
        % Update velocity
        particles.velocity(i, :) = inertia_weight * particles.velocity(i, :) + ...
            c1 * r1 .* (particles.best_position(i, :) - particles.position(i, :)) + ...
            c2 * r2 .* (global_best_position - particles.position(i, :));
        
        % Update position
        particles.position(i, :) = particles.position(i, :) + particles.velocity(i, :);
        
        % Evaluate fitness
        current_fitness = your_objective_function(particles.position(i, :), num_devices);
        
        % Update personal best
        if current_fitness < particles.best_fitness(i)
            particles.best_fitness(i) = current_fitness;
            particles.best_position(i, :) = particles.position(i, :);
        end
    end
    
    % Update global best
    [min_fitness, min_index] = min(particles.best_fitness);
    if min_fitness < global_best_fitness
        global_best_fitness = min_fitness;
        global_best_position = particles.best_position(min_index, :);
    end
    
    % Display current best fitness value for each iteration
    fprintf('Iteration %d: Best Fitness = %.4f\n', iteration, global_best_fitness);
end

% Display final result
fprintf('\nFinal Result:\n');
fprintf('Global Best Fitness = %.4f\n', global_best_fitness);
fprintf('Global Best Position = ');
disp(global_best_position);



%% Equations

% Channel Statistical Model (Gamma distribution)
% To describe the channel DC gain we use the gamma distribution
x = 1:10;
a = 13.79;
b = 0.04;
f_x = gampdf(x, a, b);

% Your objective function (modify this for your specific problem)
function value = your_objective_function(x, num_devices)
    % Compute the Rk_hut and the Pk using the functions below
    Rk_hut = 0;
    Pk = 0;
    value = Rk_hut/Pk;
end

% Average throughput of the network Rk_hat
% k is the number of the kth node
% Rk_power is the average rate of the node
% q(k) is the probability of channel access of the kth node
function Rk_hat = avThrouput (k, Rk_power, q)
    temp1 = Rk_power * q(k);
    temp2 = 1;
    for i = 1:k
        if i ~= k
            temp2 = temp2 * (1 - q(i));
        end
    end
    Rk_hat = temp1 * temp2;
end


% Average rate of the kth node Rk_power

function Rk_power = avRate (Rk, k, Pk, sigma, heta, theta)
    Xk = Xk_helper(sigma, Rk, Pk, heta, theta);
    g = gammainc(Xk, k);
    G = gamma(k);
    temp1 = g/G;
    temp2 = 1 - temp1;
    Rk_power = Rk * temp2;
end


% Xk helper function

function y = Xk_helper(sigma, Rk, Pk, heta, theta)
    numerator = 2 * pi * sigma^2 * (2*Rk - 1);
    denominator = e * (abs(heta * Pk))^2 * theta^2;
    fraction = numerator/denominator;
    y = sqrt(fraction);
end