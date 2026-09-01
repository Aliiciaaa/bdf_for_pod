function [a_val,k] = bdf1(a0, rhs, jrhs, Mr, tiempos, max_iter)
% BDF1 - Integración temporal del sistema proyectado POD utilizando el método BDF1 (Euler implícito).
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
%   a_val     - Matriz cuyos columnas contienen los coeficientes de la solución
%               aproximada en la base POD en cada instante de tiempo especificado
%               en 'tiempos'. El tamaño es [r, N], donde r es el número de modos POD
%               y N es el número de pasos de tiempo.

% Número de pasos de tiempo
N = length(tiempos);
dt = abs(tiempos(end) - tiempos(end-1));  % Paso de tiempo constante
TOL = (dt)/100;
% Inicialización de la matriz de soluciones en base POD
a_val = zeros(size(a0,1), N);

% Condición inicial (proyección sobre base POD)
a_val(:,1) = a0;

% Bucle temporal
for n = 2:N
    t_next = tiempos(n);        % Tiempo actual
    a_prev = a_val(:, n-1);     % Solución en el paso anterior
    a = a_prev;                 % Inicialización de Newton con el valor anterior

    for k = 1:max_iter
        R = Mr * (a - a_prev) / dt - rhs(t_next, a); % formulación implícita

        J = Mr / dt - jrhs(t_next, a);

        delta = -J \ R;

        % Actualización de la solución
        a = a + delta;

        if norm(delta) < TOL
            break;
        end
    end
    a_val(:, n) = a;
end
end