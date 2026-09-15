% clear all; clc;
M      = 1024; % número de pasos en tiempo
%ListNr = [18; 22; 26; 30; 34; 38; 42; 46; 50; 54];   % número de modos
% ListNr = [18; 26; 34; 42; 50];
ListNr = 18
NumNr  = length(ListNr);
erH1= zeros(NumNr,13); erL2 = erH1; 

% Descargar los autovalores del POD, y las snapshots UN (diferencias
% finitas)
load(['..\','00_compute_pod_basis\','base_pod_diff_M1024.mat']);
lambda = s.^2;
% Definimos la matriz de errores 
for i = 1:NumNr
    nr = ListNr(i);
    nombre_base = 'results_POD_';
    nombre = [nombre_base, num2str(M),'_r',num2str(nr),'.mat']; % archivo con los errores guardados

    load(nombre); % Descargamos los errores 
    % columna 1 : r (modo)
    erH1(i,1) = nr; 
    erL2(i,1) = nr;
    %
    % columna 2 : (sum cola de los autovalores)^(1/2)
    erH1(i,2) = sqrt(sum(lambda(nr+1:end))); 
    erL2(i,2) = sqrt(sum(lambda(nr+1:end)));
    %
    % columna 3 : sqrt de la media del error entre la proyección de los snapshots y los snapshots 
    % (diferencias finitas) -> coincide con la cola de los autvalores
    % ^(1/2)
    c = PhiDiff(:,1:nr)' *  kron(speye(2), Sh)* UN; 
    PrUN = PhiDiff(:,1:nr) * c; 
    EPrUN = PrUN - UN; % error proyección - snapshots
    erH1(i,3) = sqrt(mean(abs(sum(EPrUN .* (kron(speye(2), Sh) * EPrUN), 1)))); 
    erL2(i,3) = sqrt(mean(abs(sum(EPrUN .* (kron(speye(2), Mh) * EPrUN), 1))));
    % 
    % columna 4: maximo | PrYh + p - Zh| donde Zh es la aproximación FEM
    c_pry_H1 = PhiDiff(:,1:nr)'*kron(speye(2),Sh)*Yh;
    PrYhH1 = PhiDiff(:,1:nr) * c_pry_H1; 
    PrYhL2 = PhiDiff(:,1:nr) * c_pry_H1; 
    %Mc = PhiDiff(:,1:nr)'*kron(speye(2),Mh)*PhiDiff(:,1:nr);
    % coef_l2 = Mc\c_pry_L2; PrYhL22 = PhiDiff(:,1:nr)*coef_l2;
    % EPrYhL22 = PrYhL22 - Yh; 
    %erL2(i,4) = max(sqrt(abs(sum(EPrYhL22 .* (kron(speye(2), Mh) * EPrYhL22)))));
    EPrZH1 = PrYhH1 -Yh; % EPrZL2 = PrYhL2 - Yh;
    erL2PrYh = sqrt(abs(sum(EPrZH1 .* (kron(speye(2), Mh) * EPrZH1))));
    erH1PrYh = sqrt(abs(sum(EPrZH1 .* (kron(speye(2), Sh) * EPrZH1))));
    erH1(i,4) = max(erH1PrYh(5:end)); erL2(i,4) = max(erL2PrYh(5:end)); 
    %
    % columna 5 : mean | PrYh + p - Zh| donde Zh es la aproximación FEM
    erH1(i,5) = mean(erH1PrYh); erL2(i,5) = mean(erL2PrYh);
    % columna 6 : maximo |PrYh + p - Zr| donde Zr es la aproximación POD
    % con condiciones Dirichlet no homogéneas con BDF5
    erH1(i,6) = max(errH1PrYr_bdf5); 
    erL2(i,6) = max(errL2PrYr_bdf5);
    %
    % columna 7 : maximo |PrYh + p - Zr| donde Zr es la aproximación POD
    % con condiciones Dirichlet no homogéneas con ode15s
    erH1(i,7) = max(errH1PrYr_ode15s); 
    erL2(i,7) = max(errL2PrYr_ode15s);
    %
    % columna 8 - 13 : máximo |Zh - Zr| donde Zr es la aproximacón POD %
    % con condiciones Dirichlet no homogéneas con bdf q (q=1...5) y ode15s
    erH1(i,8:13) = maxerrH1h;
    erL2(i,8:13) = maxerrL2h;
end
