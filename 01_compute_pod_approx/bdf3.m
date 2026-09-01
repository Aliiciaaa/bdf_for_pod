function [a_val,k_mean] = bdf3(a0, rhs, jrhs, Mr, tiempos, max_iter)
% BDF3 - Integración temporal del sistema reducido POD utilizando el método BDF3.
%
% INPUT:
%   a0        - Vector columna con los coeficientes iniciales en la base POD
%               (proyección inicial del sistema FEM sobre la base POD).
%   rhs       - Función handle @(t, a) que define el segundo miembro (no lineal)
%               del sistema reducido en función del tiempo t y del estado a.
%   jrhs      - Función handle @(t, a) que devuelve el Jacobiano de rhs respecto
%               a 'a', necesario para el método de Newton.
%   Mr        - Matriz de masa reducida del sistema proyectado en la base POD.
%   tiempos   : tiempos en los que se calcula la solución
%   TOL       - Tolerancia para la convergencia del método de Newton-Raphson.
%   max_iter  - Número máximo de iteraciones de Newton permitidas por paso de tiempo.
%
% OUTPUT:
%   a_val     - Matriz cuyas columnas contienen los coeficientes de la solución
%               aproximada en la base POD en cada instante de tiempo especificado
%               en 'tiempos'. El tamaño es [r, N], donde r es el número de modos POD
%               y N es el número de pasos de tiempo.

% Inicialización
N = length(tiempos);                   % Número de pasos de tiempo
dt = abs(tiempos(2) - tiempos(1));    % Paso de tiempo estimado
TOL = (dt)^3/100;
t0 = tiempos(1);
t2 = tiempos(3);

a_val = zeros(size(a0,1), N);         % Arreglo para soluciones POD
a_val(:,1) = a0;                      % Asignar condición inicial

% Generar condiciones iniciales adicionales usando BDF2
% con paso refinado: dt_bdf2 = dt / q, con q entero grande
q = ceil(abs(1 / dt^(1/2)));          % Refinamiento para precisión
dt_bdf2 = dt / q;
tiempos_bdf2 = t0 : dt_bdf2 : t2;

[a_incond,~] = bdf2(a0, rhs, jrhs, Mr, tiempos_bdf2, max_iter);

% Extraer condiciones iniciales desde integración con BDF2
a_val(:,1:3) = [a_incond(:,1), a_incond(:,q+1), a_incond(:,end)];

% Ajustar tiempos para que coincidan con los puntos evaluados por BDF2
tiempos(2) = t0 + q * dt_bdf2;
tiempos(3) = t0 + 2 * q * dt_bdf2;

% Bucle temporal principal con BDF3
k_all = zeros(1,N-4);
for n = 4:N
    t_n = tiempos(n);                    % Tiempo actual
    dt = abs(tiempos(n) - tiempos(n-1)); % Paso de tiempo actual

    % Soluciones anteriores
    a_nm1 = a_val(:,n-1);  % a^{n-1}
    a_nm2 = a_val(:,n-2);  % a^{n-2}
    a_nm3 = a_val(:,n-3);  % a^{n-3}

    a = 3*a_nm1 - 3*a_nm2 + a_nm3;  % Inicialización de Newton con el último valor conocido

    for k = 1:max_iter
        % Residuo del sistema BDF3:
        % (11a^n - 18a^{n-1} + 9a^{n-2} - 2a^{n-3}) / (6dt) = rhs(t_n, a^n)
        R = Mr * ((11*a - 18*a_nm1 + 9*a_nm2 - 2*a_nm3) / (6*dt)) - rhs(t_n, a);

        % Jacobiano del sistema
        J = (11/6) * Mr / dt - jrhs(t_n, a);

        delta = -J \ R;

        % Actualización
        a = a + delta;

        if norm(delta) < TOL
            break;
        end
    end

    % Advertencia si Newton no converge
    if k == max_iter
        warning('Newton no convergió en el paso %d', n);
    end
    k_all(n) = k;
    % Guardar la solución
    a_val(:,n) = a;
end
3
k_mean = mean(k_all);
end