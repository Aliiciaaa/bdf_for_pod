% Code for computing snapshots.
clear all; close all;
nu = 0.002; alpha = 1; beta = 3; 
k = 1; % number of periods

%FEM system parameters
l = 2;   % r = 1 -> elementos lineales, r = 2 -> elementos cuadráticos
% Descargamos las condiciones iniciales, en la frontera, la triangulación y el tiempo de periodo.
load initial_conditions_nu0p002_alpha1_beta3.mat
tic;
GD = [GN;GE]; % condiciones Dirichlet
G0 = [GS;GW]; % condiciones Neumann
tic;
Tri = Tri';
% Calculamos los nodos con condiciones Dirichlet.
I_dir = unique(GD);   % Nodos con Dirichlet
nn = length(z);
I_all = (1:nn)';               % Todos los nodos
I_int = setdiff(I_all, I_dir); % Nodos interiores = nodos - Dirichlet
nint = length(I_int);
[Mh, Sh] = matrices2d(Tri,z,l);

%Coordenadas de los vértices de todos los triángulos.
x1 = z(Tri(:,1),1); x2 = z(Tri(:,2),1); x3 = z(Tri(:,3),1); 
y1 = z(Tri(:,1),2); y2 = z(Tri(:,2),2); y3 = z(Tri(:,3),2); 
% Jacobiano del cambio de variable : [xi_x xi_y; eta_x, eta_y]./J
J = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1); %determinante del jacobiano

% Definimos la condición inicial
y0 = [u(Iu);v(Iv)];

TOL = 1e-12; % tolerancia para la integración en tiempo. 
param = [nu alpha beta];
        
% qh = [int alpha * phi1; int alpha phi2, ... , int alpha phinn] = 
qh = alpha*Mh*ones(nn,1);
% Calculamos el término lineal : (con las condiciones Dirichlet incluidas)
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
% añadimos el valor en las condiciones de la frontera
Zh = zeros(2*nn,length(tiempos));
Zh(Iuv,:) = Y'; Zh(I_dir,:) = alpha; Zh(nn+I_dir,:) = beta/alpha;

computation_time_fem = toc;
nombre = 'snapshots_alpha1_beta3_nu0p002';
save(nombre, ...
     'param', 'tp', 'tiempos', 'Zh', 'Tri', 'z', 'Iu', 'Iuv', 'I_dir','Iv', 'Ahuu', 'Ahvv', 'Ahvu', 'qh', 'Mh', 'Sh','J', 'l','TOL', 'M2',...
     'S2');
% save in the same folder to compute POD basis
disp(strcat(['Snapshots saved on file',blanks(1),nombre,'.mat']))