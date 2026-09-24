% Code for getting POD approximation using finite differences and BDF-q as time integrator
close all; clear all

% Choose number of time steps for the POD basis
M = 512;
nombre_baseM= ['base_pod_diff_M',num2str(M),'.mat'];
load(['..\','00_compute_pod_basis\',nombre_baseM]);
Phi = PhiDiff;

MhnGlobal = kron(speye(2),Mh); ShnGlobal = kron(speye(2),Sh);

% Parameters set-up
nu    = 0.002;      
alpha = 1;        
beta  = 3;        
k     = 1;        % final time = k * period time 
%nrs     = [18, 26, 34, 42, 50];   % number of modes for POD method
nrs   = 18;

K = 2048/M; nn = length(z);
p = mean(Zh,2); 
Yh = Zh - p;     % FEM approximation with homogeneous boundary conditions
tiempos = 0:tp/M:tp;

tic;
S = make_tensor(PhiDiff(:,1:max(nrs)),Tri,z,I_dir,J,2,p); 
time_tensor = toc;

Ah = [Ahuu, sparse(1,1,0,nn,nn); Ahvu Ahvv];
for i = 1:length(nrs)
    nr = nrs(i);
    Sr = S(1:nr, 1:nr+1, 1:nr+1, 1:nr+1);
    get_POD_tensor;
end
