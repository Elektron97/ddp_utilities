% Copyright (C) 2019 Maitreya Venkataswamy - All Rights Reserved

% Class for quadratic terminal and running costs that inherits from
% abstract Cost class
classdef QuadraticCost < Cost
    properties
        Q_f; % terminal cost weights
        R; % control cost weights
    end
    
    methods
        % Constructor
        function obj = QuadraticCost(Q_f, R)
            obj.Q_f = Q_f;
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
        function L = L(obj, ~, u, dt)
            L = (0.5 .* u.' * obj.R * u) .* dt;
        end
        
        % Running cost derivatives
        function L_x = L_x(~, x, ~, ~)
            L_x = zeros(numel(x), 1);
        end
        function L_u = L_u(obj, ~, u, dt)
            L_u = (obj.R * u) .* dt;
        end
        function L_xx = L_xx(~, x, ~, ~)
            L_xx = zeros(numel(x));
        end
        function L_uu = L_uu(obj, ~, ~, dt)
            L_uu = obj.R .* dt;
        end
        function L_xu = L_xu(~, x, u, ~)
            L_xu = zeros(numel(x), numel(u));
        end
        function L_ux = L_ux(~, x, u, ~)
            L_ux = zeros(numel(u), numel(x));
        end
    end
end

