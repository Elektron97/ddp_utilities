% Modification by SoRoSim team ;)

%%%      DDP for GVS Hybrid-Kinematic Chain                   %%%
% This code is a readaption of the code developed by            %
% Maitreya Venkataswamy, to fit in the SoRoSim Framework.       %
% Original code: https://github.com/maitreyakv/ddp-simulation.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The idea here is to merge the solutions that we found online, mimicking
% also the algorithm described here:
% https://www.imgeorgiev.com/2023-02-01-ddp/.
% At the same time, we want to keep the structure of the optimTraj package,
% since we're working on the comparison on the two methods.

function sol = my_ddp(dyn_obj, cost_obj, x0, tf, N)
    arguments
        dyn_obj
        cost_obj
        x0
        tf
        N
    end

    %% Addpath to the Classes
    addpath("classes");

    %% Allocate variable
    % Time stamp array and timestep
    t = linspace(0.0, t_f, N);
    dt = t(2) - t(1);

    % Cost history
    J = zeros(num_iter, 1);

    % Input
    u = zeros(dyn_obj.nact, N);

    % State Space
    x = [x0, zeros(dyn_obj.ndof, N - 1)];

    % Value function and derivatives
    V = zeros(N, 1);

    % % Derivatives
    % V_x = cell(N, 1);
    % V_xx = cell(N, 1);
    
    % % State action value derivatives
    % Q_x = cell(N, 1);
    % Q_u = cell(N, 1);
    % Q_xx = cell(N, 1);
    % Q_uu = cell(N, 1);
    % Q_xu = cell(N, 1);
    % Q_ux = cell(N, 1);

    %% Nominal Input (Random Input)
    disp("initializing random control sequence...")
    
    % Random Input
    u(:, 1:end-1) = rand(dyn_obj.nact, N - 1);
    % Implicit: u(N) = 0;
    disp("generating initial trajectory...")

    %% Nominal Trajectory & Analytical Derivatives
    for i = 1:(N - 1)
        % Gradients
    end


    % % State trajectory
    % x = cell(N, 1);
    % x_new = cell(N, 1);
    % x{1} = x_0;
    % x_new{1} = x_0;
    % 
    % 
    % % Control Input trajectory
    % u = cell(N, 1);
    % 
    % % Value function and derivatives
    % V = zeros(N, 1);
    % V_x = cell(N, 1);
    % V_xx = cell(N, 1);
    % 
    % % State action value derivatives
    % Q_x = cell(N, 1);
    % Q_u = cell(N, 1);
    % Q_xx = cell(N, 1);
    % Q_uu = cell(N, 1);
    % Q_xu = cell(N, 1);
    % Q_ux = cell(N, 1);
end