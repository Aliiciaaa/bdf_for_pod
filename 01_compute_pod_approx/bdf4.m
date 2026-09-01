function [a_val,k_mean] = bdf4(a0, rhs, jrhs, Mr, tiempos, max_iter)
% BDF4 - Integración temporal del sistema proyectado POD utilizando el método BDF4.
%
% INPUTS:
%   a0        : Vector de condiciones iniciales proyectadas (coeficientes iniciales en base POD)
%   rhs       : Función del segundo miembro (RHS) evaluada como rhs(t, a)
%   jrhs      : Jacobiano del segundo miembro con respecto a a, jrhs(t, a)
%   Mr        : Matriz de masa reducida (M_r)
%   tiempos   : tiempos en los que se calcula la solución
%   TOL       : Tolerancia para la convergencia de Newton-Raphson
%   max_iter  : Número máximo de iteraciones de Newton por paso de tiempo
%
% OUTPUT:
%   a_val     : Matriz de coeficientes reducidos en cada paso de tiempo

    N = length(tiempos);                      % Número total de pasos de tiempo
    dt = abs(tiempos(2) - tiempos(1));        % Paso de tiempo asumido uniforme
    TOL = (dt)^4/100;
    t0 = tiempos(1);
    t3 = tiempos(4);                          % Tiempo hasta el cual se hace integración previa con BDF3

    a_val = zeros(size(a0,1), N);             % Inicializa la matriz de soluciones

    % BDF4 necesita 4 pasos previos -> se obtiene con pasos más pequeños de BDF3

    % Se calcula el número de subpasos para BDF3 para igualar integración hasta t3
    q = ceil(abs(1 / dt^(1/3)));              % Escoge q para que (dt_bdf3)^3 ~ dt
    dt_bdf3 = dt / q;
    tiempos_bdf3 = t0:dt_bdf3:t3;

    [a_incond, ~] = bdf3(a0, rhs, jrhs, Mr, tiempos_bdf3, max_iter);
    a_val(:,1:4) = [a_incond(:,1), a_incond(:,q+1), a_incond(:,2*q+1), a_incond(:,end)];

    % Corregir los tiempos intermedios usados en BDF4 con los tiempos reales
    tiempos(2) = t0 + q * dt_bdf3;
    tiempos(3) = t0 + 2 * q * dt_bdf3;
    tiempos(4) = t0 + 3 * q * dt_bdf3;
    
    k_all = zeros(1,N-5);
    for n = 5:N
        t_n = tiempos(n);
        dt = abs(tiempos(n) - tiempos(n-1));

        a_nm1 = a_val(:,n-1);
        a_nm2 = a_val(:,n-2);
        a_nm3 = a_val(:,n-3);
        a_nm4 = a_val(:,n-4);

        % Inicialización para Newton-Raphson : extrapolación 
        a = 4*a_nm4 -6*a_nm3 + 4*a_nm2 - a_nm1;

        for k = 1:max_iter
            R = Mr * ( (25*a - 48*a_nm1 + 36*a_nm2 - 16*a_nm3 + 3*a_nm4) / (12*dt) ) - rhs(t_n, a);
            J = (25/12) * Mr / dt - jrhs(t_n, a);

            delta = -J \ R;
            a = a + delta;

            if norm(delta) < TOL
                break;
            end
        end

        if k == max_iter
            warning('Newton no convergió en el paso de tiempo %d (t = %.5f)', n, t_n);
        end
        k_all(n) = k;
        a_val(:,n) = a;
    end
4
k_mean = mean(k_all);
end