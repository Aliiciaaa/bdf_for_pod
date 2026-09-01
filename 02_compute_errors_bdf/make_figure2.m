% This code creates Figure 2 in the article.
clear all;
close all;

nombre_base2 = 'results_POD_';
num_metodos = 5;
Mss = [64,128,256,512,1024];
nM = length(Mss);
nr = 18;
errL2 = zeros(num_metodos, nM);
errH1 = zeros(num_metodos, nM);

for j = 1:nM
    M = Mss(j);
    nombre2 = [nombre_base2, num2str(M),'_r',num2str(nr),'.mat'];
    filePath = fullfile('..','01_compute_pod_approx', nombre2);
    load(filePath)

    maxerrL2 = [max(abs(errL2CombinedBdf1(2:end))), ...
                max(abs(errL2CombinedBdf2(3:end))), ...
                max(abs(errL2CombinedBdf3(4:end))), ...
                max(abs(errL2CombinedBdf4(5:end))), ...
                max(abs(errL2CombinedBdf5(6:end)))];

    maxerrH1 = [max(abs(errH1CombinedBdf1(2:end))), ...
                max(abs(errH1CombinedBdf2(3:end))), ...
                max(abs(errH1CombinedBdf3(4:end))), ...
                max(abs(errH1CombinedBdf4(5:end))), ...
                max(abs(errH1CombinedBdf5(6:end)))];

    errL2(:, j) = maxerrL2;
    errH1(:, j) = maxerrH1;
end

colores = lines(num_metodos);

%% === L2 ERROR ===
figL2 = figure('Units','centimeters','Position',[5 5 16 10]); 
for n = 1:num_metodos
    grid off;
    loglog(Mss, errL2(n,:), ...
        'Color', colores(n,:), ...
        'LineStyle','-', ...
        'Marker','x', ...
        'DisplayName',['BDF' num2str(n)], ...
        'LineWidth',1.5);

    p = polyfit(log10(Mss(3:end)), log10(errL2(n,3:end)), 1);
    slope = p(1);

    x_text = Mss(end);
    y_text = errL2(n,end);

    label_str = ['$\mathrm{slope} = ' num2str(slope,'%.2f') '$'];
    text(x_text*1.05, y_text*1.15, label_str, ...
        'FontSize',11, ...
        'Interpreter','latex', ...
        'Color',colores(n,:), ...
        'FontWeight','bold');

    hold on
end
hold off

xlabel('$M$','Interpreter','latex','FontSize',14)
title('Errors in $L^2$','Interpreter','latex','FontSize',14)
legend('show','Location','southwest', ...
       'FontSize',11,'Interpreter','latex')

set(gca,'FontSize',12,'TickLabelInterpreter','latex')

xlim([0 2200])
ylim([10^(-7) 10])
xticks(Mss)
grid on

exportgraphics(figL2,'maxerrL2.pdf','ContentType','vector')


%% === H1 ERROR ===
figH1 = figure('Units','centimeters','Position',[5 5 16 10]); 
for n = 1:num_metodos
    grid off;
    loglog(Mss, errH1(n,:), ...
        'Color', colores(n,:), ...
        'LineStyle','-', ...
        'Marker','x', ...
        'DisplayName',['BDF' num2str(n)], ...
        'LineWidth',1.5);

    p = polyfit(log10(Mss(3:end)), log10(errL2(n,3:end)), 1);
    slope = p(1);

    x_text = Mss(end);
    y_text = errH1(n,end);

    label_str = ['$\mathrm{slope} = ' num2str(slope,'%.2f') '$'];
    text(x_text*1.05, y_text*1.15, label_str, ...
        'FontSize',11, ...
        'Interpreter','latex', ...
        'Color',colores(n,:), ...
        'FontWeight','bold');

    hold on
end
hold off

xlabel('$M$','Interpreter','latex','FontSize',14)
title('Errors in $H^1_0$','Interpreter','latex','FontSize',14)
legend('show','Location','southwest', ...
       'FontSize',11,'Interpreter','latex')

set(gca,'FontSize',12,'TickLabelInterpreter','latex')

xlim([0 2200])
ylim([10^(-6) 10])
xticks(Mss)
grid on

exportgraphics(figH1,'maxerrH1.pdf','ContentType','vector')