% Code for computing finite differences in time.
clear all; close all;

%Ms = [32,64,128,256,512]; % number of time steps
Ms = 1024;
beta = 3; alpha = 1;
load('snapshots_alpha1_beta3_nu0p002.mat'); % load snapshots for the corresponding beta and alpha

for i = 1 : length(Ms)

    M = Ms(i);

    % Extraemos U y V desde UV
    K = 2048/M;
    Zhk = Zh(:,1:K:end); p = mean(Zhk,2);
    
    % Yh tiene condiciones Dirichlet homogéneas.
    % anohter option is to take p as
    % nn = size(Zh,1)/2;
    % p1 = ones(nn,1)*alpha; p2 = ones(nn,1)*beta/alpha;
    % p = [p1;p2];

    Yh = Zhk - p;
    
    dt = tp/M;
    if size(Yh, 2) < 2
        error('No hay suficientes columnas para calcular la derivada temporal.');
    end
    Yt = diff(Yh, [], 2) / dt;
    
    % Construcción de la matriz UN
    
    UN = Yt; % we do not consider the initial term
    
    name = ['snapshots_diff_M',num2str(M),'.mat'];
    save(name);
    
    disp(['Resultados guardados en el archivo ', name, '.mat']);
end