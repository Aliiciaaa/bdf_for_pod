function [a_val,k_mean] = bdf4(a0, rhs, jrhs, Mr, tiempos, max_iter)
        % BDF4 - Temporal integration of the projected POD system using the BDF4 method.
        %
        % INPUTS:
        %   a0        : Vector of projected initial conditions (initial coefficients in the POD basis)
        %   rhs       : Right-hand side (RHS) function evaluated as rhs(t, a)
        %   jrhs      : Jacobian of the right-hand side with respect to a, jrhs(t, a)
        %   Mr        : Reduced mass matrix (M_r)
        %   tiempos   : Time steps at which the solution is computed
        %   TOL       : Tolerance for the Newton-Raphson convergence
        %   max_iter  : Maximum number of Newton iterations per time step
        %
        % OUTPUT:
        %   a_val     : Matrix of reduced coefficients at each time step

        N = length(tiempos);                      
dt = abs(tiempos(2) - tiempos(1));        
TOL = (dt)^4 / 100;
t0 = tiempos(1);
t3 = tiempos(4); % Time up to which the preliminary integration with BDF3 is performed

a_val = zeros(size(a0, 1), N);            

% BDF4 requires 4 previous steps -> obtained using smaller steps of BDF3

% Calculate the number of substeps for BDF3 to match the integration up to t3
q = ceil(abs(1 / dt^(1/3))); % Choose q so that (dt_bdf3)^3 ~ dt
dt_bdf3 = dt / q;
tiempos_bdf3 = t0:dt_bdf3:t3;

% Correct the intermediate times used in BDF4 with the actual times
[a_incond, ~] = bdf3(a0, rhs, jrhs, Mr, tiempos_bdf3, max_iter);
a_val(:, 1:4) = [a_incond(:, 1), a_incond(:, q+1), a_incond(:, 2*q+1), a_incond(:, end)];

tiempos(2) = t0 + q * dt_bdf3;
tiempos(3) = t0 + 2 * q * dt_bdf3;
tiempos(4) = t0 + 3 * q * dt_bdf3;

k_all = zeros(1, N); % Se inicializa con tamaño N completo para evitar problemas de indexación mas adelante

for n = 5:N
    t_n = tiempos(n);
    dt = abs(tiempos(n) - tiempos(n-1));

    a_nm1 = a_val(:, n-1);
    a_nm2 = a_val(:, n-2);
    a_nm3 = a_val(:, n-3);
    a_nm4 = a_val(:, n-4);

    % Predictor inicial para el método de Newton
    a = 4*a_nm4 - 6*a_nm3 + 4*a_nm2 - a_nm1;

    for k = 1:max_iter
        R = Mr * ((25*a - 48*a_nm1 + 36*a_nm2 - 16*a_nm3 + 3*a_nm4) / (12*dt)) - rhs(t_n, a);
        J = (25/12) * Mr / dt - jrhs(t_n, a);

        delta = -J \ R;
        a = a + delta;

        if norm(delta) < TOL
            break;
        end
    end

    if k == max_iter
        warning('Newton did not converge at step %d (BDF4) [t = %.5f]', n, t_n);
    end
    k_all(n) = k;
    a_val(:, n) = a;
end
4
k_mean = mean(k_all(5:end));
end
