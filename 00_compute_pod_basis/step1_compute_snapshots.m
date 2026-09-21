% Code for computing snapshots.
clear all; close all;
nu = 0.002; alpha = 1; beta = 3; 
k = 1; % number of periods

%FEM system parameters
l = 2;   % l = 1 -> linear elements, l = 2 -> quadratic elements
% Load initial conditions, boundary points, triangulation and period time.
load initial_conditions_nu0p002_alpha1_beta3.mat

GD = [GN;GE]; % condiciones Dirichlet
G0 = [GS;GW]; % condiciones Neumann
Tri = Tri';

I_dir = unique(GD);            % Nodes with Dirichlet conditions
nn = length(z);                % Number of nodes
I_all = (1:nn)';               % All nodes
I_int = setdiff(I_all, I_dir); % Interior nodes = All nodes - Dirichlet nodes
nint = length(I_int);

[Mh, Sh] = matrices2d(Tri,z,l);

x1 = z(Tri(:,1),1); x2 = z(Tri(:,2),1); x3 = z(Tri(:,3),1); 
y1 = z(Tri(:,1),2); y2 = z(Tri(:,2),2); y3 = z(Tri(:,3),2); 
% Change of variable for the Jacobian : [xi_x xi_y; eta_x, eta_y]./J
J = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1); 

y0 = [u(Iu);v(Iv)];

TOL = 1e-12; % % time integrator tolerance
param = [nu alpha beta];
        
% qh = [int alpha * phi1; int alpha phi2, ... , int alpha phinn] = 
qh = alpha*Mh*ones(nn,1);
Ahuu = (-nu*Sh - (beta+1)*Mh);
Ahvv = -nu*Sh;
Ahvu = beta*Mh;      
Ah = [Ahuu, sparse(1,1,0,nn,nn); Ahvu, Ahvv];

rhs = @(t,y) rhsq2d(t, y, Tri, z, I_dir, Ahuu, Ahvv, Ahvu, qh, J, l, param);
jrhs = @(t,y) jrhsq2d(t, y, Tri, z, I_dir,Ahuu, Ahvv, Ahvu, J, l, param);

tiempos = 0:tp/2048:tp;
Iuv = [Iu; nn + Iv];

M2 = kron(speye(2),Mh);
M2_int = kron(speye(2),Mh(Iu,Iu));
S2 = kron(speye(2),Sh);
opes = odeset('AbsTol', TOL/1000,'RelTol', TOL,'Jacobian', jrhs, 'Mass', M2_int, 'Stats','on');
[tiempos_ode,Y] = ode15s(rhs, tiempos , y0 , opes);

Zh = zeros(2*nn,length(tiempos));  % add boundary values
Zh(Iuv,:) = Y'; Zh(I_dir,:) = alpha; Zh(nn+I_dir,:) = beta/alpha;

nombre = 'snapshots_alpha1_beta3_nu0p002';
save(nombre, ...
     'param', 'tp', 'tiempos', 'Zh', 'Tri', 'z', 'Iu', 'Iuv', 'I_dir','Iv', 'Ahuu', 'Ahvv', 'Ahvu', 'qh', 'Mh', 'Sh','J', 'l','TOL', 'M2',...
     'S2');
disp(strcat(['Snapshots saved on file',blanks(1),nombre,'.mat']))
