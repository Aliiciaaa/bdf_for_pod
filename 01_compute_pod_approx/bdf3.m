function [a_val,k_mean] = bdf3(a0, rhs, jrhs, Mr, tiempos, max_iter)
        % BDF3 - Temporal integration of the reduced POD system using the BDF3 method.
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
TOL = (dt)^3/100;
t0 = tiempos(1);
t2 = tiempos(3);

a_val = zeros(size(a0,1), N);         
a_val(:,1) = a0;                     

% Generate additional initial conditions using BDF2
% with a refined time step: dt_bdf2 = dt / q, where q is a large integer
q = ceil(abs(1 / dt^(1/2)));         
dt_bdf2 = dt / q;
tiempos_bdf2 = t0 : dt_bdf2 : t2;

[a_incond,~] = bdf2(a0, rhs, jrhs, Mr, tiempos_bdf2, max_iter);

% Extract initial conditions from the BDF2 integration
a_val(:,1:3) = [a_incond(:,1), a_incond(:,q+1), a_incond(:,end)];

% Adjust times to match the points evaluated by BDF2
tiempos(2) = t0 + q * dt_bdf2;
tiempos(3) = t0 + 2 * q * dt_bdf2;

k_all = zeros(1,N-4);
for n = 4:N
    t_n = tiempos(n);                    
    dt = abs(tiempos(n) - tiempos(n-1)); 

    a_nm1 = a_val(:,n-1);  % a^{n-1}
    a_nm2 = a_val(:,n-2);  % a^{n-2}
    a_nm3 = a_val(:,n-3);  % a^{n-3}

    a = 3*a_nm1 - 3*a_nm2 + a_nm3; 

    for k = 1:max_iter
        % (11a^n - 18a^{n-1} + 9a^{n-2} - 2a^{n-3}) / (6dt) = rhs(t_n, a^n)
        R = Mr * ((11*a - 18*a_nm1 + 9*a_nm2 - 2*a_nm3) / (6*dt)) - rhs(t_n, a);

        J = (11/6) * Mr / dt - jrhs(t_n, a);

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
    a_val(:,n) = a;
end
3
k_mean = mean(k_all);
end
