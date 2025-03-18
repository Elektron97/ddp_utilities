% Copyright (C) 2019 Maitreya Venkataswamy - All Rights Reserved

% Abstract class for system dynamics
classdef (Abstract) Dynamics
    methods
        %% Equations of Motion
        dxdt = F(obj, t, x, u)
        
        %% Analytical Derivatives
        [fx, fu] = analytical_derivatives(obj, t, x, xdot, u);
    end
end