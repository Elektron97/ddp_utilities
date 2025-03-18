% Copyright (C) 2019 Maitreya Venkataswamy - All Rights Reserved

% Abstract class for system dynamics
classdef (Abstract) Dynamics
    methods
        %% Equations of Motion
        x_dot = dynamics(obj, t, x, u)
        
        %% Analytical Derivatives
        [fx, fu] = analytical_derivatives(obj, t, x, xdot, u)
    end
end