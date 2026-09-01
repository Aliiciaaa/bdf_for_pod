function [a_val,k_mean] = bdf2(a0, rhs, jrhs, Mr, tiempos, max_iter)
% BDF2 - Integración temporal del sistema reducido POD utilizando el método BDF2.
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
N = length(tiempos);                     % Número total de pasos de tiempo
dt = abs(tiempos(2) - tiempos(1));       % Paso de tiempo inicial
TOL = (dt)^2/100;
t0 = tiempos(1); 
t1 = tiempos(2);

a_val = zeros(size(a0,1), N);            % Soluciones POD

% Paso inicial usando BDF1 para obtener dos primeros valores
% Se usa paso de tiempo refinado: dt_bdf1 = dt^2
dt_bdf1 = dt;
tiempos_bdf1 = t0:dt_bdf1:t1;

[a_incond,~] = bdf1(a0, rhs, jrhs, Mr, tiempos_bdf1, max_iter);
a_val(:,1:2) = [a_incond(:,1), a_incond(:,end)];

% Ajuste del segundo valor de tiempo para consistencia con bdf1
tiempos(2) = tiempos_bdf1(end);

% Integración en tiempo con BDF2
k_all = zeros(1,N-3);
for n = 3:N
    t_n = tiempos(n);         % Tiempo actual
    dt = abs(tiempos(n) - tiempos(n-1)); % Paso de tiempo actual

    % Valores anteriores
    a_nm1 = a_val(:,n-1);     % a^{n-1}
    a_nm2 = a_val(:,n-2);     % a^{n-2}

    a = 2*a_nm1 - a_nm2;                % Inicialización de Newton

    for k = 1:max_iter
        % Evaluación del residuo de la ecuación implícita
        % (3a^n - 4a^{n-1} + a^{n-2}) / (2dt) = rhs(t_n, a^n)
        R = Mr * ((3*a - 4*a_nm1 + a_nm2) / (2*dt)) - rhs(t_n, a);

        % Evaluación del Jacobiano
        J = (3/2) * Mr / dt - jrhs(t_n, a);

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
    % Guardar solución actual
    a_val(:,n) = a;
end
2
k_mean = mean(k_all);
end