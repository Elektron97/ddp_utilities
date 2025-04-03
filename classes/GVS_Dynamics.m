%%% Class for a generic SoRoSim object, that inherits the abstract Dynamics
%%% class to fit the DDP problem.

classdef GVS_Dynamics < Dynamics
    properties
        robot_linkage       % SoRoSim Object
        ndof                % n. of DoFs
        nx                  % dimension of the state space
        nact
    end

    methods
        %% Constructor
        function obj = GVS_Dynamics(robot_linkage)
            % Init
            obj.robot_linkage = robot_linkage;
            obj.ndof = obj.robot_linkage.ndof;
            obj.nx = 2*obj.ndof;
            obj.nact = obj.robot_linkage.nact;
        end

        %% Dynamics
        function x_dot = dynamics(obj, t, x, u)
            % dynamicsSolver function
            [y,~,~,~] = obj.robot_linkage.dynamicsSolver(t, x, u);
            
            % dxdt = [q_dot; q_2dot]
            x_dot = [x(obj.ndof + 1:end); y(1:obj.ndof)];
        end

        %% Analytical Derivatives (First Order)
        function [fx, fu] = analytical_derivatives(obj, t, x, xdot, u, dt)
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

            %% Implement discretization
            fx = eye(obj.nx) + fx*dt;
            fu = fu*dt;
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
                [fx1, ~] = obj.analytical_derivatives(0, x_tilde1, obj.dynamics(0, x_tilde1, u), u, dt);
                [fx2, ~] = obj.analytical_derivatives(0, x_tilde2, obj.dynamics(0, x_tilde2, u), u, dt);

                % Finite Difference
                fxx(:, :, i) = 0.5*(fx2 - fx1)/options.eps;
            end

            % Differentiation w.r.t. u
            for i = 1:obj.nact
                % Perturbating x
                u_tilde1 = u - options.eps*u_bases(:, i);
                u_tilde2 = u + options.eps*u_bases(:, i);

                % Compute Analytical Derivatives
                [fx1, fu1] = obj.analytical_derivatives(0, x, obj.dynamics(0, x, u_tilde1), u_tilde1, dt);
                [fx2, fu2] = obj.analytical_derivatives(0, x, obj.dynamics(0, x, u_tilde2), u_tilde2, dt);

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
