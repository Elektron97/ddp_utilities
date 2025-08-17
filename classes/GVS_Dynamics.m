%%% Class for a generic SoRoSim object, that inherits the abstract Dynamics
%%% class to fit the DDP problem.

classdef GVS_Dynamics < Dynamics
    properties
        robot_linkage       % SoRoSim Object
        ndof                % n. of DoFs
        nx                  % dimension of the state space
        nact
        fsolve_opt          % Store fsolve options for Implicit Time Integration
    end

    methods
        %% Constructor
        function obj = GVS_Dynamics(robot_linkage)
            % Init
            obj.robot_linkage = robot_linkage;
            obj.ndof = obj.robot_linkage.ndof;
            obj.nx = 2*obj.ndof;
            obj.nact = obj.robot_linkage.nact;

            % Init fsolve options
            obj.fsolve_opt = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg', ...
                                            'Display', 'iter', 'SpecifyObjectiveGradient', true);
        end

        %% Dynamics
        function x_dot = dynamics(obj, t, x, u)
            % dynamicsSolver function
            [y,~,~,~] = obj.robot_linkage.dynamicsSolver(t, x, u);
            
            % dxdt = [q_dot; q_2dot]
            x_dot = [x(obj.ndof + 1:end); y(1:obj.ndof)];
        end

        %% Discretization of Dynamics
        function x_new = discretize_dynamics(obj, t, x, u, h, options)
            arguments
                obj
                t
                x
                u
                h
                options.method = "ode1"
                options.specify_jacobian = false
            end

            %% Discretization
            switch(options.method)
                case "ode1"
                    % Compute Continous Dynamics
                    x_dot = obj.dynamics(t, x, u);

                    x_new = x + h*x_dot;
                case "ode4"
                    % Compute Continous Dynamics
                    x_dot = obj.dynamics(t, x, u);

                    k0 = x_dot;
                    k1 = obj.dynamics(t + 0.5*h, x + 0.5*h*k0, u);
                    k2 = obj.dynamics(t + 0.5*h, x + 0.5*h*k1, u);
                    k3 = obj.dynamics(t + h, x + h*k2, u);
                    
                    % Update State
                    x_new = x + (h/6)*(k0 + 2*k1 + 2*k2 + k3);
                case "ode45"
                    ODEFUN = @(t, xk) obj.dynamics(t, xk, u);
                    [~, YOUT] = ode45(ODEFUN, [t, t + h], x);
                    x_new = YOUT(end, :)';
                case "ode15s"
                    ODEFUN = @(t, xk) obj.dynamics(t, xk, u);

                    % Increase Speed
                    if options.specify_jacobian
                        odeopts = odeset('Jacobian', @(t, xk) obj.analytical_derivatives(t, xk, obj.dynamics(t, xk, u), u));
                        [~, YOUT] = ode15s(ODEFUN, [t, t + h], x, odeopts);
                    else
                        [~, YOUT] = ode15s(ODEFUN, [t, t + h], x);
                    end

                    x_new = YOUT(end, :)';
                
                %% Implicit Time Integration
                case "ode1i"
                    % Create Function handle of the Residual
                    ODEFUN = @(x_next) computeResidualAndJacobian(t, x_next, x, u, h, obj, "method", options.method);

                    % fsolve: Initial Guess x_k (continuity)
                    x_new = fsolve(ODEFUN, x, obj.fsolve_opt);  

                case "trapz"
                    % Create Function handle of the Residual
                    ODEFUN = @(x_next) computeResidualAndJacobian(t, x_next, x, u, h, obj, "method", options.method);

                    % fsolve: Initial Guess x_k (continuity)
                    x_new = fsolve(ODEFUN, x, obj.fsolve_opt);

                otherwise
                    error("Unsupported Integration Method.");
            end
        end

        %% Analytical Derivatives (First Order)
        function [fx, fu] = analytical_derivatives(obj, t, x, xdot, u)
            q    = x(1:obj.ndof);
            qd   = x(obj.ndof +1:2*obj.ndof);
            qdd  = xdot(obj.ndof +1:2*obj.ndof);
            
            %% DAE Jacobian
            dae_index = 1;
            % lambda is required in case of closed-loop joint
            lambda = [];
            [~,~,~,dID_dq,dID_dqd,dID_dqdd,~,dtau_dq,dtau_dqd,dtau_du,~,~,~,~,~] = obj.robot_linkage.DAEJacobians(t, q, qd, qdd, u, lambda, dae_index);
            
            % Computing Derivatives
            B = dtau_du;
            M = dID_dqdd;
            dFDdq  = M\(dtau_dq-dID_dq);
            dFDdqd = M\(dtau_dqd-dID_dqd);

            %% Gradients
            fx = [zeros(obj.ndof, obj.ndof), eye(obj.ndof); dFDdq, dFDdqd];
            fu = [zeros(obj.ndof, obj.nact); M\B];
        end

        %% Apply Discretization
        function [fx, fu] = discretized_gradients(obj, t, x, u, h, options)
            arguments
                obj
                t
                x
                u
                h
                options.method = "ode1"
            end

            %% Gradients of the Discretized System
            switch(options.method)
                case "ode1"
                    % Compute Analytical Derivatives (Continous System)
                    xdot = obj.dynamics(t, x, u);
                    [fx_continous, fu_continous] = obj.analytical_derivatives(t, x, xdot, u);

                    fx = eye(obj.nx) + fx_continous*h;
                    fu = fu_continous*h;
                case "ode4"
                    % Compute Analytical Derivatives (Continous System)
                    xdot = obj.dynamics(t, x, u);
                    [fx_continous, fu_continous] = obj.analytical_derivatives(t, x, xdot, u);

                    k0 = xdot;  % xdot = f(x, u) | Continous
                    k1 = obj.dynamics(t + 0.5*h, x + 0.5*h*k0, u);
                    k2 = obj.dynamics(t + 0.5*h, x + 0.5*h*k1, u);
                    k3 = obj.dynamics(t + h, x + h*k2, u);

                    % Compute Analytical Derivatives of i_th k
                    % Init
                    k0x = fx_continous;
                    k0u = fu_continous;

                    % k-ith
                    [k1x, k1u] = obj.analytical_derivatives(t + 0.5*h, x + 0.5*h*k0, k1, u);
                    [k2x, k2u] = obj.analytical_derivatives(t + h/2, x + 0.5*h*k1, k2, u);
                    [k3x, k3u] = obj.analytical_derivatives(t + h, x + h*k2, k3, u);

                    % Chain Rule Correction (Input) | diff_u(x_k + 0.5*h*k0(u))
                    k1u = k1u + 0.5*h*k1x*k0u;
                    k2u = k2u + 0.5*h*k2x*k1u;
                    k3u = k3u + h*k3x*k2u;

                    % Chain Rule Correction (State) | diff_x(x_k + 0.5*h*k0) = 
                    % = eye() + 0.5*h*k0x.
                    k1x = k1x*(eye(obj.nx) + 0.5*h*k0x);
                    k2x = k2x*(eye(obj.nx) + 0.5*h*k1x);
                    k3x = k3x*(eye(obj.nx) + h*k2x);

                    % Rewrite Gradients
                    fx = eye(obj.nx) + (h/6)*(k0x + 2*k1x + 2*k2x + k3x);
                    fu = (h/6)*(k0u + 2*k1u + 2*k2u + k3u);

                %% Implicit Time Integration
                case "ode1i"
                    % Compute Analytical Jacobian
                    x_new = obj.discretize_dynamics(t, x, u, h, "method", options.method);
                    [~, gx_next, gx, gu] = computeResidualAndJacobian(t, x_new, x, u, h, obj, "method", options.method);
                    fx = - gx_next\gx;
                    fu = - gx_next\gu;
                
                case "trapz"
                    % Compute Analytical Jacobian
                    x_new = obj.discretize_dynamics(t, x, u, h, "method", options.method);
                    [~, gx_next, gx, gu] = computeResidualAndJacobian(t, x_new, x, u, h, obj, "method", options.method);
                    fx = - gx_next\gx;
                    fu = - gx_next\gu;

                otherwise
                    error("Unsupported Discretization Method.");
            end

        end

        %% Numerical Hessian
        function [fxx, fxu, fuu] = numerical_hessian(obj, x, u, dt, options)
            arguments
                obj
                x
                u
                dt
                options.eps = 1e-6
            end
            
            % Init Hessians
            fxx = zeros(obj.nx, obj.nx, obj.nx);
            fxu = zeros(obj.nx, obj.nx, obj.nact);
            fuu = zeros(obj.nx, obj.nact, obj.nact);

            % Storing Canonical Base
            x_bases = eye(obj.nx);
            u_bases = eye(obj.nact);

            % Differentiation w.r.t. x
            for i = 1:obj.nx
                % Perturbating x
                x_tilde1 = x - options.eps*x_bases(:, i);
                x_tilde2 = x + options.eps*x_bases(:, i);

                % Compute Analytical Derivatives
                [fx1, ~] = obj.discretized_gradients(0, x_tilde1, obj.dynamics(0, x_tilde1, u), u, dt);
                [fx2, ~] = obj.discretized_gradients(0, x_tilde2, obj.dynamics(0, x_tilde2, u), u, dt);

                % Finite Difference
                fxx(:, :, i) = 0.5*(fx2 - fx1)/options.eps;
            end

            % Differentiation w.r.t. u
            for i = 1:obj.nact
                % Perturbating x
                u_tilde1 = u - options.eps*u_bases(:, i);
                u_tilde2 = u + options.eps*u_bases(:, i);

                % Compute Analytical Derivatives
                [fx1, fu1] = obj.discretized_gradients(0, x, obj.dynamics(0, x, u_tilde1), u_tilde1, dt);
                [fx2, fu2] = obj.discretized_gradients(0, x, obj.dynamics(0, x, u_tilde2), u_tilde2, dt);

                % Finite Difference
                fxu(:, :, i) = 0.5*(fx2 - fx1)/options.eps;
                fuu(:, :, i) = 0.5*(fu2 - fu1)/options.eps;
            end

            %% Symmetrize
            fxx = 0.5*(fxx + permute(fxx, [1 3 2]));
            fuu = 0.5*(fuu + permute(fuu, [1 3 2]));
        end

    end
end

function [F, Jx_next, Jx, Ju] = computeResidualAndJacobian(t, x_next, x_current, u, h, dyn_obj, options)
    arguments
        t
        x_next
        x_current
        u
        h
        dyn_obj
        options.method = "ode1i"     % obj is the object containing your dynamics and derivatives methods
    end

    % f(x_next, u)
    x_next_dot = dyn_obj.dynamics(t, x_next, u);

    % Compute only if fsolve requires Jacobian
    if nargout > 1
        [fx_next, fu_next] = dyn_obj.analytical_derivatives(t, x_next, x_next_dot, u);
    end

    switch(options.method)
        case "ode1i"
            % 1. Compute Residual
            F = x_next - x_current - h * x_next_dot;

            if nargout > 1
                % 2. Compute the Jacobian of the residual
                n = numel(x_current); % Get the size of the state vector
                I = eye(n);
                
                Jx_next = I - h * fx_next;
                Jx = - I;
                Ju = - h * fu_next;
            end

        case "trapz"
            % f(x_current, u)
            x_current_dot = dyn_obj.dynamics(t, x_current, u);

            % 1. Compute Residual
            F = x_next - x_current - 0.5*h*(x_next_dot + x_current_dot);

            if nargout > 1
                [fx_current, fu_current] = dyn_obj.analytical_derivatives(t, x_current, x_current_dot, u);

                % 2. Compute the Jacobian of the residual
                n = numel(x_current); % Get the size of the state vector
                I = eye(n);
                
                Jx_next = I - 0.5*h*fx_next;
                Jx = - I - 0.5*h*fx_current;
                Ju = - 0.5*h*(fu_current + fu_next);
            end

        otherwise
            error("Unsupported Implicit Integration Method.")
    end
end

%% Newton - Raphson Method
function [x, exitflag, output] = newton_raphson(res_and_jac_func, x0, options)
    % NEWTON_RAPHSON Solves a system of nonlinear equations F(x) = 0.
    %
    %   This function finds a root of a system of nonlinear equations using the
    %   Newton-Raphson method. It is optimized to solve the linear system
    %   J*dx = -F at each iteration rather than computing a matrix inverse.
    %
    %   INPUTS:
    %       res_and_jac_func - Function handle. Must take a vector x and return
    %                          [F, J], where F is the residual vector F(x) and
    %                          J is the Jacobian matrix dF/dx.
    %       x0               - Column vector of the initial guess.
    %       options          - (Optional) Struct with solver options:
    %                          .Tol (1e-8)  : Tolerance for the norm of the residual.
    %                          .MaxIter (50): Maximum number of iterations.
    %
    %   OUTPUTS:
    %       x        - The solution vector (the root).
    %       exitflag - An integer indicating the reason the algorithm terminated:
    %                  1: Converged to the specified tolerance.
    %                  0: Reached the maximum number of iterations without converging.
    %       output   - A struct with details about the solution process:
    %                  .iterations: Number of iterations performed.
    %                  .final_norm_F: The norm of the residual at the solution.
    
    arguments
        res_and_jac_func    function_handle
        x0                  (:,1) double
        options.Tol         (1,1) double = 1e-6
        options.MaxIter     (1,1) double = 1400
    end
    
    % Initialize variables
    x = x0;
    iter = 0;
    
    % Start the Newton-Raphson iterations
    for iter = 1:options.MaxIter
        % Get the residual and the Jacobian from the user-provided function
        [F, J] = res_and_jac_func(x);
        
        % Check for convergence
        norm_F = norm(F);
        if norm_F < options.Tol
            exitflag = 1; % Success
            output.iterations = iter - 1;
            output.final_norm_F = norm_F;
            return;
        end
        
        % --- Core Newton Step ---
        % Solve the linear system J*delta_x = -F for the update step delta_x.
        % This is far more numerically stable and computationally efficient
        % than calculating delta_x = inv(J) * -F.        
        % Update the solution guess
        x = x - J \ F;
    end
    
    % If the loop completes, max iterations were reached without convergence
    exitflag = 0; % Failure
    output.iterations = options.MaxIter;
    output.final_norm_F = norm(F);
end
