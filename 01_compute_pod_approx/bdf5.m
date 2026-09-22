function [a_val,k_mean] = bdf5(a0, rhs, jrhs, Mr, tiempos, max_iter)
        % BDF5 - Temporal integration of the projected POD system using the BDF5 method 
        %        with Newton-Raphson iteration at each time step.
        % INPUTS:
        %   a0        : Column vector with the initial condition (initial projection onto the POD basis)
        %   rhs       : Function evaluating the right-hand side of the reduced system: rhs(t, a)
        %   jrhs      : Function evaluating the Jacobian of the right-hand side with respect to 'a': jrhs(t, a)
        %   Mr        : Mass matrix of the reduced system (e.g., Phi' * Mh * Phi)
        %   tiempos   : Time steps at which the solution is computed
        %   TOL       : Tolerance for the Newton-Raphson convergence criterion
        %   max_iter  : Maximum number of Newton-Raphson iterations per time step
        %
        % OUTPUT:
        %   a_val     : Matrix containing the coefficients of the projected solution (POD basis) at each time instant
        %
        % NOTE:
        %   The first five time steps are solved using the BDF4 method, with a finer
        %   temporal subdiscretization, to generate the required initial conditions for BDF5.

N = length(tiempos);
dt = abs(tiempos(2) - tiempos(1));
TOL = (dt)^5/100;
t0 = tiempos(1); 
t4 = tiempos(5);

a_val = zeros(size(a0,1), N);

% To obtain the first five values required by BDF5, BDF4 is used with a refined time step
q = ceil(abs(1/dt^(1/4)));
dt_bdf4 = dt / q;
tiempos_bdf4 = t0:dt_bdf4:t4;

% Solve with BDF4 and extract the required points
[a_incond, ~] = bdf4(a0, rhs, jrhs, Mr, tiempos_bdf4, max_iter);
a_val(:,1:5) = [a_incond(:,1), a_incond(:,q+1), a_incond(:,2*q+1), a_incond(:,3*q+1), a_incond(:,end)];

% Adjust the times to match the interpolated values from BDF4
tiempos(2) = tiempos_bdf4(q+1); 
tiempos(3) = tiempos_bdf4(2*q+1);
tiempos(4) = tiempos_bdf4(3*q+1);
tiempos(5) = tiempos_bdf4(4*q+1);

% === BDF5 from step 6 ===
k_all = zeros(1,N-6);
for n = 6:N
    t_n = tiempos(n);

    a_nm1 = a_val(:,n-1);
    a_nm2 = a_val(:,n-2);
    a_nm3 = a_val(:,n-3);
    a_nm4 = a_val(:,n-4);
    a_nm5 = a_val(:,n-5);

    a = 5*a_nm5 -10*a_nm4 +10*a_nm3 -5*a_nm2 + a_nm1;

    for k = 1:max_iter
        R = Mr * ( (137*a - 300*a_nm1 + 300*a_nm2 - 200*a_nm3 + 75*a_nm4 - 12*a_nm5) / (60*dt) ) - rhs(t_n, a);
        
        J = (137/60) * Mr / dt - jrhs(t_n, a);
        
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
    a_val(:,n) = a;
end
5
k_mean = mean(k_all);
end
