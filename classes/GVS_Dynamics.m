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

    end
end
