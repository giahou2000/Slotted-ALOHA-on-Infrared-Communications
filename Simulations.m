%% Performance of Slotted ALOHA on infrared WBANs
% This simulation tries to describe the performance of Slotted ALOHA
% when applied as a MAC protocol for a WBAN system that consists of
% 6 nodes (5 sensor nodes and 1 coordinator node). For the communication medium it is assumed that we use the
% infrared part of the electromagnetic spectrum.

% Authors:
% Giachoudis Christos, christos.giachoudis@fresnel.fr
% Vasilis Papanikolaou,
% Konstantinos Rallis,



%% PSO - Particle Swarm Optimization
% A brute force / exploration method for the optimization of energy
% efficiency in Wireless Body Area Networks



%% !___ PSO Parameters definition and initialization ___!

% Limits on the variable parameters of the problem

clc; clear;

qLow = 0; % probability
PLow = 0; % (mw)
RLow = 250; % (kbps)

qHigh = 1; % probability
PHigh = 30; % (mW)
RHigh = 2500; % (kbps)

% The values blo and bup represent the lower and upper boundaries of the search-space respectively.

blo = [qLow, PLow, RLow]; % minimum limit for the tested variables
bup = [qHigh, PHigh, RHigh]; % maximum limit for the tested variables

% Velocity limits

velocity_low = -abs(bup - blo);
velocity_high = abs(bup - blo);
num_devices = 5; % Change this to the desired number of devices (K)
num_variables = 3;
num_particles = 30; % number of particles that will search for the best position
max_iterations = 100; % iterations until an acceptable convergence
phi_p = 1.5; % cognitive parameter
phi_g = 2; % social parameter
w = 0.7; % inertia weight
global_best_fitness = 0;
s = 0.5;
heta = 0.6;
theta = 0.2;

% ______________________________________________

% Initialize particles
% Pre-allocate empty arrays of structs
particles = struct.empty(num_particles, 0);  % 0 indicates no pre-defined fields
global_best = struct.empty(num_devices, 0);

% Access and modify individual structs
for i = 1:num_particles
    for j = 1:num_devices
      particles(i).node(j).position(1) = qLow + (qHigh - qLow) * rand; % Random initial probability
      particles(i).node(j).position(2) = PLow + (PHigh - PLow) * rand; % Random initial power
      particles(i).node(j).position(3) = RLow + (RHigh - RLow) * rand; % Random initial transmission rate
      particles(i).node(j).velocity = velocity_low + (velocity_high - velocity_low) * rand; % Random initial velocities
      particles(i).node(j).best_position = particles(i).node(j).position; % Best known positions
      particles(i).best_fitness = 0; % Best known fitness values
    end
end

% Evaluate starting fitness for each particle
for i = 1:num_particles
    particles(i).best_fitness = fitness(particles, num_devices, i, s, heta, theta);
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





%% Main PSO loop
% ______________________________________________

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
            if 5>3 % implement the limits of the variables in if statement
                particles(i).node(j).position = temp;
            else
                particles(i).node(j).position = old_position; % instead of old position put the edge
            end
        end
        
        % Evaluate fitness
        current_fitness = fitness(particles, num_devices, i, s, heta, theta);
        
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

% Display final result
fprintf('\nFinal Result:\n');
fprintf('Global Best Fitness = %.4f\n', global_best_fitness);
fprintf('Global Best Position = ');
for i = 1:num_devices
    disp(global_best(i).position);
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
function value = fitness(particles, nodes_num, which_particle, s, heta, theta)
    % value = rand; % just for testing the functionality of the rest of the
    % algorithm
    % Compute the Rk_hut and the Pk using the functions below
    sum = 0;
    for i = 1:nodes_num
        Rk_power = avRate (particles(which_particle).node(i).position(3), i, particles(which_particle).node(i).position(2), s, heta, theta);
        Rk_hat = avThrouput (i, Rk_power, particles(which_particle).node(i).position(1));
        Pk = particles(which_particle).node(i).position(2);
        sum = sum + (Rk_hat/Pk);
    end
    value = sum;
end

% Average throughput of the network Rk_hat
% k is the number of the kth node
% Rk_power is the average rate of the node
% q(k) is the probability of channel access of the kth node
function Rk_hat = avThrouput (k, Rk_power, q)
    temp1 = Rk_power * q;
    temp2 = 1;
    for i = 1:k
        if i ~= k
            temp2 = temp2 * (1 - q);
        end
    end
    Rk_hat = temp1 * temp2;
end


% Average rate of the kth node Rk_power

function Rk_power = avRate (Rk, k, Pk, s, heta, theta)
    Xk = Xk_helper(s, Rk, Pk, heta, theta);
    g = gammainc(Xk, k);
    G = gamma(k);
    temp1 = g/G;
    temp2 = 1 - temp1;
    Rk_power = Rk * temp2;
end


% Xk helper function

function y = Xk_helper(s, Rk, Pk, heta, theta)
    numerator = 2 * pi * (s^2) * ((2 * Rk) - 1);
    denominator = exp(1) * (abs(heta * Pk))^2 * theta^2;
    fraction = numerator/denominator;
    if fraction < 0
        y = 0;
    else
        y = sqrt(fraction);
    end
end