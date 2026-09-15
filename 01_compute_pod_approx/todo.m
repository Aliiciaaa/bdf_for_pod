% Aproximación POD usando diferencias finitas con las bdf como integrador temporal.
close all; clear all

M = 64;
nombre_baseM= ['base_pod_diff_M',num2str(M),'.mat'];
load(['..\','00_compute_pod_basis\',nombre_baseM]);
Phi = PhiDiff;
MhnGlobal = kron(speye(2),Mh); ShnGlobal = kron(speye(2),Sh);

% Parámetros.
nu    = 0.002;      
alpha = 1;        
beta  = 3;        
k     = 1;        % número de periodos (tiempo final = k * tiempo del periodo)
nrs     = [18, 26, 34, 42, 50];       % modos para el método POD

K = 2048/M; nn = length(z);
p = mean(Zh,2); 
Yh = Zh - p;
tiempos = 0:tp/M:tp;
tic;
S = make_tensor(PhiDiff(:,1:max(nrs)),Tri,z,I_dir,J,2,p); 
time_tensor = toc;

Ah = [Ahuu, sparse(1,1,0,nn,nn); Ahvu Ahvv];
for i = 1:length(nrs)
    nr = nrs(i);
    Sr = S(1:nr, 1:nr+1, 1:nr+1, 1:nr+1);
    get_POD_tensor;
    disp('... done.');
end