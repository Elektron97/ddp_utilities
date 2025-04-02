% Copyright (C) 2019 Maitreya Venkataswamy - All Rights Reserved 

% Modification by SoRoSim team
%%%      DDP for GVS Hybrid-Kinematic Chain                   %%%
% This code is a readaption of the code developed by            %
% Maitreya Venkataswamy, to fit in the SoRoSim Framework.       %
% Original code: https://github.com/maitreyakv/ddp-simulation.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% We improved the code using https://github.com/imgeorgiev/ddp %%%
% https://github.com/Alescontrela/DDP_Examples.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Differential Dynamic Programming Algorithm Function

% Performs Differential Dynamic Programming Algorithm to produce a
% locally optimal control sequence.
%
% Inputs
% 
% x_0      : initial state vector
% x_star   : target state vector
% t_f      : time horizon
% N        : number of discretizations of time
% dyn      : instance of Dynamics class
% cost     : instance of Cost class
% u_max    : maximum magnitude of control, used for clamping
% num_iter : number of DDP iterations
% alpha    : DDP learning parameter for line searching
%
% Outputs
% 
% sol : structure with solution components
%
function sol = ddp_backtracking(x0, x_goal, t_f, N, dyn, cost, num_iter, options)
    arguments
        x0
        x_goal
        t_f
        N
        dyn
        cost
        num_iter
        options.show_online = true
        options.cost_tolerance = 1e-5
        options.backtracking = true
        options.bt_it = 10
        options.bt_decay = 0.5
        options.reg = 1e-6
        options.u0 = (-1.0 + 2*rand(dyn.nact, N))
    end

    %% Init
    % Time stamp array and timestep
    t = linspace(0.0, t_f, N);
    dt = t(2) - t(1);

    % Cost history
    J = zeros(num_iter, 1);

    % State trajectory
    x = [x0, zeros(dyn.nx, N - 1)];

    % Store x_dot for Analytical Derivatives
    x_dot = zeros(dyn.nx, N);
    % the last will be 0

    % Input
    u = zeros(dyn.nact, N);

    % Value function and derivatives
    % V = zeros(1, N);
    V_x = zeros(dyn.nx, N);
    V_xx = zeros(dyn.nx, dyn.nx, N);

    % State action value derivatives
    Q_x = zeros(dyn.nx, N);
    Q_u = zeros(dyn.nact, N);
    Q_xx = zeros(dyn.nx, dyn.nx, N);
    Q_uu = zeros(dyn.nact, dyn.nact, N);
    Q_xu = zeros(dyn.nx, dyn.nact, N);
    Q_ux = zeros(dyn.nact, dyn.nx, N);

    % Alphas value
    alphas = options.bt_decay.^(0:(options.bt_it - 1));

    %% Nominal Trajectory
    % Random (for now)
    u = options.u0;

    % Propagation with the Dynamics
    for k = 1:N-1
        x_dot(:, k) = dyn.dynamics(t(k), x(:, k), u(:, k));
        x(:, k + 1) = x(:, k) + dt*x_dot(:, k);
    end

    %% Show Online
    if options.show_online
        f1 = figure;
        plot(nan, nan);
    end

    %% Main Loop
    for i = 1:num_iter
        disp("Iteration " + num2str(i));

        %% Backwards pass
        [V_x, V_xx, Q_x, Q_u, Q_xx, Q_uu, Q_xu, Q_ux] = backward_pass(dyn, cost, x_goal, t, x, u, x_dot, N, dt, options.reg);

        %% Forward Pass
        if options.backtracking && i > 1
            % Backtracking
            for j = 1:options.bt_it
                % Forward Pass
                [x_star, u_star, x_dot_star] = forward_pass(dyn, t, x, u, Q_u, Q_ux, Q_uu, N, dt, alphas(j));
    
                % Compute the new cost
                J_try = compute_totalCost(cost, x_star, u_star, x_goal, N, dt);
    
                % Check
                if(J_try < J(i - 1))
                    break;
                end
            end

            % Update
            x = x_star;
            u = u_star;
            x_dot = x_dot_star;
            J(i) = J_try;
        else
            %% Forward Pass
            [x, u, x_dot] = forward_pass(dyn, t, x, u, Q_u, Q_ux, Q_uu, N, dt, alphas(1));

            %% Compute total cost
            J(i) = compute_totalCost(cost, x, u, x_goal, N, dt);
        end

        %% Stopping Conditions
        % Error for NaN
        if(anynan(x) || anynan(u))
            disp("Error: NaN in the solutions. Stopping...");
            % Assemble solution
            sol = assemble_solution(x, u, t, J, Q_u, Q_uu, Q_ux, 1);
            break;
        end

        %% Show Total Cost
        if i > 1
            if options.show_online
                figure(f1)
                plot(J(1:i), 'LineWidth', 2.0, 'Marker', 'o')
                grid on
                xlabel("N° of iterations")
                ylabel("Total Cost J")
                xlim([1, i])
                drawnow
            end
    
            % Cost
            if(abs(J(i) - J(i - 1)) < options.cost_tolerance)
                disp("Stop DDP since the difference of the costs are lower then the tolerance.");
                % Assemble solution
                sol = assemble_solution(x, u, t, J, Q_u, Q_uu, Q_ux, 0);
                break;
            end
        end
    end  

    %% Assemble and return solution structure
    disp("finished DDP, assembling results for post-processing...\n");

    % Assemble solution
    sol = assemble_solution(x, u, t, J, Q_u, Q_uu, Q_ux, 0);
end

%% Utilities
function J = compute_totalCost(cost, x, u, x_goal, N, dt)
    J = 0;
    for k = 1:N-1
        J = J + cost.L(x(:, k), u(:, k), dt);
    end
    J = J + cost.phi(x(:, N), x_goal);
end

function [V_x, V_xx, Q_x, Q_u, Q_xx, Q_uu, Q_xu, Q_ux] = backward_pass(dyn, cost, x_goal, t, x, u, x_dot, N, dt, reg)
    % Init
    % Value function and derivatives
    % V = zeros(1, N);
    V_x = zeros(dyn.nx, N);
    V_xx = zeros(dyn.nx, dyn.nx, N);

    % State action value derivatives
    Q_x = zeros(dyn.nx, N);
    Q_u = zeros(dyn.nact, N);
    Q_xx = zeros(dyn.nx, dyn.nx, N);
    Q_uu = zeros(dyn.nact, dyn.nact, N);
    Q_xu = zeros(dyn.nx, dyn.nact, N);
    Q_ux = zeros(dyn.nact, dyn.nx, N);

    % Compute terminal value function and derivatives
    % V(N) = cost.phi(x(:, N), x_goal);
    V_x(:, N) = cost.phi_x(x(:, N), x_goal);
    V_xx(:, :, N) = cost.phi_xx(x(:, N), x_goal);

    % Perform backwards pass
    for k = N-1:-1:1
        % Compute Analytical Derivatives
        [fx, fu] = dyn.analytical_derivatives(t(k), x(:, k), x_dot(:, k), u(:, k), dt);

        Q_x(:, k) = cost.L_x(x(:, k), u(:, k), dt) + fx.' * V_x(:, k + 1);
        Q_u(:, k) = cost.L_u(x(:, k), u(:, k), dt) + fu.' * V_x(:, k + 1);
        Q_xx(:, :, k) = cost.L_xx(x(:, k), u(:, k), dt) + fx.' * V_xx(:, :, k + 1) * fx;
        Q_uu(:, :, k) = cost.L_uu(x(:, k), u(:, k), dt) + fu.' * V_xx(:, :, k + 1) * fu;
        Q_xu(:, :, k) = cost.L_xu(x(:, k), u(:, k), dt) + fx.' * V_xx(:, :, k + 1) * fu;
        Q_ux(:, :, k) = cost.L_ux(x(:, k), u(:, k), dt) + fu.' * V_xx(:, :, k + 1) * fx;

        % Compute the value function derivatives
        % Regularization
        Q_uu(:, :, k) = Q_uu(:, :, k) + reg*eye(size(Q_uu(:, :, k)));
    
        % Update Gradients
        V_x(:, k) = Q_x(:, k) - Q_xu(:, :, k) * (Q_uu(:, :, k) \ Q_u(:, k));
        V_xx(:, :, k) = Q_xx(:, :, k) - Q_xu(:, :, k) * (Q_uu(:, :, k) \ Q_ux(:, :, k));
    end
end

function [x_star, u_star, x_dot_star] = forward_pass(dyn, t, x, u, Q_u, Q_ux, Q_uu, N, dt, alpha)
    % Init        
    x_star = zeros(size(x));
    x_dot_star = zeros(size(x));
    x_star(:, 1) = x(:, 1);
    u_star = zeros(size(u));

    % Backtracking
    for k = 1:N-1
        error = x_star(:, k) - x(:, k);
        
        % % WrapToP Correction
        % error(1:2) = wrapToPi(error(1:2));

        % Update control
        u_star(:, k) = u(:, k) - (Q_uu(:, :, k)) \ (alpha.*Q_u(:, k) + Q_ux(:, :, k)*error);

        % Compute next state in trajectory with new control
        % Store x_dot
        x_dot_star(:, k) = dyn.dynamics(t(k), x_star(:, k), u_star(:, k));
        x_star(:, k + 1) = x_star(:, k) + x_dot_star(:, k) .* dt;
    end
end

%% Helper Functions for DDP Algorithm
% Assembles and returns solution structure
%
% Inputs
%
% x               : locally optimal state trajectory
% u               : locally optimal control sequence
% t               : discretized time stamps
% J               : iteration history of cost function
% Q_u, Q_uu, Q_ux : derivatives of state action value function
% error           : zero if no error, nonzero else
%
% Outputs
%
% sol : solution structure
function sol = assemble_solution(x, u, t, J, Q_u, Q_uu, Q_ux, error)

    % Solution structure
    sol = struct;
    sol.error = error;
    sol.x = x;
    sol.u = u;
    sol.t = t;
    sol.dt = t(2) - t(1);
    sol.J = J;
    sol.Q_u = Q_u;
    sol.Q_uu = Q_uu;
    sol.Q_ux = Q_ux;
    return
end
