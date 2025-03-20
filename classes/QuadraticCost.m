% Copyright (C) 2019 Maitreya Venkataswamy - All Rights Reserved

% Class for quadratic terminal and running costs that inherits from
% abstract Cost class
classdef QuadraticCost < Cost
    properties
        Q_f;    % Terminal Cost Weights
        Q;      % State Cost Weights
        R;      % Control Cost Weights
    end
    
    methods
        % Constructor
        function obj = QuadraticCost(Q_f, Q, R)
            obj.Q_f = Q_f;
            obj.Q = Q;
            obj.R = R;
        end
        
        % Terminal state cost
        function phi = phi(obj, x_f, x_star)
            error_f = x_f - x_star;
            % WrapToPi Correction on theta1 and theta2
            error_f(1:2) = wrapToPi(error_f(1:2));

            phi = 0.5 .* (error_f).' * obj.Q_f * (error_f);
        end
        
        % Terminal state cost derivatives
        function phi_x = phi_x(obj, x_f, x_star)
            error_f = x_f - x_star;
            % WrapToPi Correction on theta1 and theta2
            error_f(1:2) = wrapToPi(error_f(1:2));

            phi_x = obj.Q_f * error_f;
        end
        function phi_xx = phi_xx(obj, ~, ~)
            phi_xx = obj.Q_f;
        end
        
        % Running cost
        function L = L(obj, x, u, dt)
            L = (0.5 .* x.' * obj.Q * x) .* dt + (0.5 .* u.' * obj.R * u) .* dt;
        end
        
        % Running cost derivatives
        function L_x = L_x(obj, x, u, dt)
            L_x = obj.Q*x*dt;
        end
        function L_u = L_u(obj, x, u, dt)
            L_u = (obj.R * u) .* dt;
        end
        function L_xx = L_xx(obj, x, u, dt)
            L_xx = obj.Q*dt;
        end
        function L_uu = L_uu(obj, x, u, dt)
            L_uu = obj.R .* dt;
        end
        function L_xu = L_xu(obj, x, u, dt)
            L_xu = zeros(numel(x), numel(u));
        end
        function L_ux = L_ux(obj, x, u, dt)
            L_ux = zeros(numel(u), numel(x));
        end
    end
end

