%%% Class for a generic SoRoSim object, that inherits the abstract Dynamics
%%% class to fit the DDP problem.

classdef GVS_Dynamics < Dynamics
    properties
        robot_linkage       % SoRoSim Object
        ndof                % n. of DoFs
        nact
    end

    methods
        %% Constructor
        function obj = GVS_Dynamics(robot_linkage)
            % Init
            obj.robot_linkage = robot_linkage;
            obj.ndof = obj.robot_linkage.ndof;
            obj.nact = obj.robot_linkage.nact;

            % Libraries
            addpath("functions")
        end

        %% Dynamics
        function x_dot = dynamics(obj, x, u)
            % For now: time-invariant implementation
            t = 0;

            % dynamicsSolver function
            [y,~,~,~] = obj.robot_linkage.dynamicsSolver(t, x, u);
            
            % dxdt = [q_dot; q_2dot]
            x_dot = [x(obj.ndof + 1:end); y(1:obj.ndof)];
        end

        %% Analytical Derivatives (First Order)
        function [fx, fu] = analytical_derivatives(obj, t, x, xdot, u)
            
        end

    end
end
