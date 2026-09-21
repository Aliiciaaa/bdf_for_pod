% Code for computing finite differences in time.
clear all; close all;

Ms = [32,64,128,256,512,1024]; % number of time steps
beta = 3; alpha = 1;
load('snapshots_alpha1_beta3_nu0p002.mat'); % load snapshots for the corresponding beta and alpha

for i = 1 : length(Ms)

    M = Ms(i);
    K = 2048/M;
    Zhk = Zh(:,1:K:end); p = mean(Zhk,2);
    
    % Yh - FEM solution with homogenenous Dirichlet conditions
    % option 2. take p as
    % nn = size(Zh,1)/2;
    % p1 = ones(nn,1)*alpha; p2 = ones(nn,1)*beta/alpha;
    % p = [p1;p2];
    
    Yh = Zhk - p; % option 1. take p as the mean
    
    dt = tp/M;
    if size(Yh, 2) < 2
        error('Unsufficient columns for computing finite differences.');
    end
    Yt = diff(Yh, [], 2) / dt;
    
    UN = Yt; % Note we do not consider the initial term
    
    name = ['snapshots_diff_M',num2str(M),'.mat'];
    save(name);
    
    disp(['Finite differences saved in ', name, '.mat']);
end
