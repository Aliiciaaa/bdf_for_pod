function [a_val,k_mean] = bdf5(a0, rhs, jrhs, Mr, tiempos, max_iter)
% BDF5 - Integración temporal del sistema proyectado POD utilizando el método BDF5 
%        utilizando iteración de Newton-Raphson en cada paso en tiempo.
% INPUTS:
%   a0        : Vector columna con la condición inicial (proyección inicial en base POD)
%   rhs       : Función que evalúa el lado derecho del sistema reducido: rhs(t, a)
%   jrhs      : Función que evalúa la Jacobiana del lado derecho respecto a 'a': jrhs(t, a)
%   Mr        : Matriz de masa del sistema reducido (por ejemplo, Phi' * Mh * Phi)
%   tiempos   : tiempos en los que se calcula la solución
%   TOL       : Tolerancia para el criterio de convergencia de Newton-Raphson
%   max_iter  : Número máximo de iteraciones de Newton-Raphson por paso de tiempo
%
% OUTPUT:
%   a_val     : Matriz que contiene los coeficientes de la solución proyectada (base POD) en cada instante de tiempo
%
% NOTA:
%   Los primeros cinco pasos temporales se resuelven utilizando el método BDF4,
%   con una subdiscretización temporal más fina, para generar las condiciones iniciales necesarias para BDF5.

% Número total de pasos de tiempo
N = length(tiempos);

% Paso de tiempo
dt = abs(tiempos(2) - tiempos(1));
TOL = (dt)^5/100;
t0 = tiempos(1); 
t4 = tiempos(5);

% Inicialización de la matriz solución
a_val = zeros(size(a0,1), N);

% Para obtener los cinco primeros valores requeridos por BDF5, se usa BDF4 con paso de tiempo refinado
q = ceil(abs(1/dt^(1/4))); % Refinamiento temporal
dt_bdf4 = dt / q;
tiempos_bdf4 = t0:dt_bdf4:t4;

% Resolver con BDF4 y extraer los puntos necesarios
[a_incond, ~] = bdf4(a0, rhs, jrhs, Mr, tiempos_bdf4, max_iter);
a_val(:,1:5) = [a_incond(:,1), a_incond(:,q+1), a_incond(:,2*q+1), a_incond(:,3*q+1), a_incond(:,end)];

% Ajuste de los tiempos para que coincidan con los valores interpolados de BDF4
tiempos(2) = tiempos_bdf4(q+1); 
tiempos(3) = tiempos_bdf4(2*q+1);
tiempos(4) = tiempos_bdf4(3*q+1);
tiempos(5) = tiempos_bdf4(4*q+1);

% === BDF5 a partir del paso 6 ===
k_all = zeros(1,N-6);
for n = 6:N
    t_n = tiempos(n);

    a_nm1 = a_val(:,n-1);
    a_nm2 = a_val(:,n-2);
    a_nm3 = a_val(:,n-3);
    a_nm4 = a_val(:,n-4);
    a_nm5 = a_val(:,n-5);

    % Inicialización para Newton
    a = 5*a_nm5 -10*a_nm4 +10*a_nm3 -5*a_nm2 + a_nm1;

    % Iteración de Newton-Raphson
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
        warning('Newton no convergió en el paso %d', n);
    end
    k_all(n) = k;
    a_val(:,n) = a;
end
5
k_mean = mean(k_all);
end