function [a_val,k_mean] = bdf2(a0, rhs, jrhs, Mr, tiempos, max_iter)
        % BDF2 - Temporal integration of the reduced POD system using the BDF2 method.
        %
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
dt = abs(tiempos(2) - tiempos(1));       
TOL = (dt)^2/100;
t0 = tiempos(1); 
t1 = tiempos(2);

a_val = zeros(size(a0,1), N);            

% Initial step using BDF1 to obtain the first two values
% A refined time step is used: dt_bdf1 = dt^2

dt_bdf1 = dt;
tiempos_bdf1 = t0:dt_bdf1:t1;

[a_incond,~] = bdf1(a0, rhs, jrhs, Mr, tiempos_bdf1, max_iter);
a_val(:,1:2) = [a_incond(:,1), a_incond(:,end)];

% Adjustment of the second time value for consistency with BDF1
tiempos(2) = tiempos_bdf1(end);

k_all = zeros(1,N-3);
for n = 3:N
    t_n = tiempos(n);         
    dt = abs(tiempos(n) - tiempos(n-1)); 

    a_nm1 = a_val(:,n-1);     % a^{n-1}
    a_nm2 = a_val(:,n-2);     % a^{n-2}

    a = 2*a_nm1 - a_nm2;                

    for k = 1:max_iter
        % (3a^n - 4a^{n-1} + a^{n-2}) / (2dt) = rhs(t_n, a^n)
        R = Mr * ((3*a - 4*a_nm1 + a_nm2) / (2*dt)) - rhs(t_n, a);

        J = (3/2) * Mr / dt - jrhs(t_n, a);

        delta = -J \ R;
        a = a + delta;

        if norm(delta) < TOL
            break;
        end
    end

    if k == max_iter
        warning('Newton did not converge at step %d', n);
    end
    k_all(n) = k;
    % Guardar solución actual
    a_val(:,n) = a;
end
2
k_mean = mean(k_all);
end
