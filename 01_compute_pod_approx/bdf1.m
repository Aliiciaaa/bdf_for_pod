function [a_val,k] = bdf1(a0, rhs, jrhs, Mr, tiempos, max_iter)
        % BDF1 - Temporal integration of the projected POD system using the BDF1 method (implicit Euler).
        % INPUT:
        %   a0        - Column vector with the initial coefficients in the POD basis
        %               (initial projection of the FEM system onto the POD basis).
        %   rhs       - Function handle @(t, a) that defines the right-hand side (non-linear)
        %               of the reduced system as a function of time t and state a.
        %   jrhs      - Function handle @(t, a) that returns the Jacobian of rhs with respect
        %               to 'a', required for the Newton method.
        %   Mr        - Reduced mass matrix of the projected system in the POD basis.
        %   tiempos   - Time steps at which the solution is computed.
        %   TOL       - Tolerance for the convergence of the Newton-Raphson method.
        %   max_iter  - Maximum number of Newton iterations allowed per time step.
        %
        % OUTPUT:
        %   a_val     - Matrix whose columns contain the coefficients of the approximate
        %               solution in the POD basis at each time instant specified
        %               in 'tiempos'. The size is [r, N], where r is the number of POD modes
        %               and N is the number of time steps.

N = length(tiempos);
dt = abs(tiempos(end) - tiempos(end-1));  % Paso de tiempo constante
TOL = (dt)/100;

a_val = zeros(size(a0,1), N);

% Initial condition 
a_val(:,1) = a0;

for n = 2:N
    t_next = tiempos(n);        % at time t_{n}
    a_prev = a_val(:, n-1);     % solution at t_{n-1}
    a = a_prev;                 % start Newton with solution at t_{n-1}
    for k = 1:max_iter
        R = Mr * (a - a_prev) / dt - rhs(t_next, a); 

        J = Mr / dt - jrhs(t_next, a);

        delta = -J \ R;
        a = a + delta;

        if norm(delta) < TOL
            break;
        end
    end
    a_val(:, n) = a;
end
end
