font = 18; fna = 18; % FontSize 

% Load POD basis
load base_pod_diff;

TOL = 1e-12;
max_iter = 20;

nombre_base = 'snaps_beta';
nombre_base2 = 'results_POD_';

if abs(beta - fix(beta)) == 0
    str_beta = ['0', num2str(beta)];
else
    str_beta = ['0', num2str(100 * beta)];
end

nombre = [nombre_base, str_beta]; % para descargar los snapshots
nombre2 = [nombre_base2, num2str(M) '.mat']; % para guardar los resultados

%% Cargar snapshots
load(nombre); % Contiene UV y tiempos, Ahuu, Ahvv, Ahvu, fh
K = 2048/M;
Uh = UV(:, 1:K:length(tiempos));
Vh = UV(:, length(tiempos)+1 :K: end);
Zh = [Uh; Vh]; p = mean(Zh,2);
Yh = Zh - p;
tiempos = 0:tp/M:tp;

Ah = [Ahuu, sparse(1,1,0,nn,nn); Ahvu Ahvv];
%% Aproximación POD
param = [nu, alpha, beta];
% with ode as time integrator
[U_ode,V_ode, a_ode] = calculate_pod_approximation(Yh(:,1),Tri,z, I_dir, J, PhiDiff(:,1:nr), Mh, Sh, Ah, fh, p, tiempos, r, TOL, param);
[U_bdf1, V_bdf1, a_bdf1] = bdf1(y0, PhiDiff, Mh, Sh, Lh, r, tiempos, TOL, max_iter, param); % la condición inicial para bdf1 es la proyección en H1
[U_bdf2, V_bdf2, a_bdf2] = bdf2(a_ode(:,1:2), PhiDiff, Mh, Lh, r, tiempos, TOL, max_iter, param); % la condición inicial para bdf2 son las dos primeras aproximaciones POD
[U_bdf3, V_bdf3, a_bdf3] = bdf3(a_ode(:,1:3), PhiDiff, Mh, Lh, r, tiempos, TOL, max_iter, param);
[U_bdf4, V_bdf4, a_bdf4] = bdf4(a_ode(:,1:4), PhiDiff, Mh, Lh, r, tiempos, TOL, max_iter, param);
[U_bdf5, V_bdf5, a_bdf5] = bdf5(a_ode(:,1:5), PhiDiff, Mh, Lh, r, tiempos, TOL, max_iter, param);

% errors in L2 and H1 between ode and bdf
[~, ~, ~, ~, errL2CombinedBdf1, errH1CombinedBdf1] = calculate_errors(U_bdf1, V_bdf1, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2Combined, errH1Combined] = calculate_errors(U_ode, V_ode, Uh, Vh, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf2, errH1CombinedBdf2] = calculate_errors(U_bdf2, V_bdf2, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf3, errH1CombinedBdf3] = calculate_errors(U_bdf3, V_bdf3, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf4, errH1CombinedBdf4] = calculate_errors(U_bdf4, V_bdf4, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf5, errH1CombinedBdf5] = calculate_errors(U_bdf5, V_bdf5, U_ode, V_ode, Mh, Sh);

maxerrL2 = [max(abs(errL2CombinedBdf1)),max(abs(errL2CombinedBdf2)),max(abs(errL2CombinedBdf3)),max(abs(errL2CombinedBdf4)),max(abs(errL2CombinedBdf5))];
maxerrH1 = [max(abs(errH1CombinedBdf1)),max(abs(errH1CombinedBdf2)),max(abs(errH1CombinedBdf3)),max(abs(errH1CombinedBdf4)),max(abs(errH1CombinedBdf5))];
% save results
save(nombre2);
disp(['Results saved in file ', nombre2]);

%% code for representing all errors wrt time
% in L2 norm
figure;

errL2CombinedBdf1(errL2CombinedBdf1 == 0) = eps; % en la escala logarítmica los ceros no se representan, cambiamos cero por eps para que aparezcan
semilogy(tiempos, errL2CombinedBdf1, 'b--', 'LineWidth', 1.5, 'DisplayName','bdf1');
hold on;
errL2CombinedBdf2(errL2CombinedBdf2 == 0) = eps;
semilogy(tiempos, errL2CombinedBdf2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'bdf2');
errL2CombinedBdf3(errL2CombinedBdf3 == 0) = eps;
semilogy(tiempos, errL2CombinedBdf3, 'r--', 'LineWidth', 1.5, 'DisplayName', 'bdf3');
errL2CombinedBdf4(errL2CombinedBdf4 == 0) = eps;
semilogy(tiempos, errL2CombinedBdf4, 'm--', 'LineWidth', 1.5, 'DisplayName', 'bdf4');
errL2CombinedBdf5(errL2CombinedBdf5 == 0) = eps;
semilogy(tiempos, errL2CombinedBdf5, 'y--', 'LineWidth', 1.5, 'DisplayName', 'bdf5');

xlabel('t [s]', 'FontSize', fna, 'Interpreter', 'latex');
ylabel('error ODE & BDF', 'FontSize', fna, 'Interpreter', 'latex');
title(['Errors in $L^2$ (M = ', num2str(M), ', r = ', num2str(r), ') with $\beta =$ ', num2str(beta, '%.2f')], 'FontSize', font, 'Interpreter', 'latex');
legend('show','Interpreter', 'LaTeX', 'FontSize', 8, 'Location', 'southeast'); 
grid on;

set(gca, 'FontSize', fna);
saveas(gcf, 'errorL2.png');

% in H1 norm
figure;

errH1CombinedBdf1(errH1CombinedBdf1 == 0) = eps;
semilogy(tiempos, errH1CombinedBdf1, 'b--', 'LineWidth', 1.5, 'DisplayName','bdf1'); 
hold on;
errH1CombinedBdf2(errH1CombinedBdf2 == 0) = eps;
semilogy(tiempos, errH1CombinedBdf2, 'g--', 'LineWidth', 1.5, 'DisplayName','bdf2'); 
errH1CombinedBdf3(errH1CombinedBdf3 == 0) = eps;
semilogy(tiempos, errH1CombinedBdf3, 'r--', 'LineWidth', 1.5, 'DisplayName','bdf3'); 
errH1CombinedBdf4(errH1CombinedBdf4 == 0) = eps;
semilogy(tiempos, errH1CombinedBdf4, 'm--', 'LineWidth', 1.5, 'DisplayName', 'bdf4');
errH1CombinedBdf5(errH1CombinedBdf5 == 0) = eps;
semilogy(tiempos, errH1CombinedBdf5, 'y--', 'LineWidth', 1.5, 'DisplayName', 'bdf5');

xlabel('t [s]', 'FontSize', fna, 'Interpreter', 'latex');
ylabel('error ODE & BDF', 'FontSize', fna, 'Interpreter', 'latex');
title(['Errors in $H^1$ (M = ', num2str(M), ', r = ', num2str(r), ') with $\beta =$ ', num2str(beta, '%.2f')], 'FontSize', font, 'Interpreter', 'latex');
legend('show','Interpreter', 'LaTeX', 'FontSize', 8, 'Location', 'southeast'); 
grid on;
set(gca, 'FontSize', fna);
saveas(gcf, 'errorH1.png');