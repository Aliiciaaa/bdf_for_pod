
font = 18; fna = 18; % FontSize 

max_iter = 20; % number maximum of iterations for bdf schemes

name_results = 'results_POD_tensor_';
full_name_results = [name_results, num2str(M),'_r', num2str(nr),'.mat']; % to save results

a0  = Phi(:,1:nr)' * ShnGlobal * Yh(:,1); % initial condition
Mr  = Phi(:,1:nr)' * kron(speye(2), Mh) * Phi(:,1:nr); % Mass matrix POD system 
Ar  = Phi(:,1:nr)' * Ah * Phi(:,1:nr); % Stiffness matrix POD system 
Arp = Phi(:,1:nr)' * Ah * p;
qr  = Phi(:,1:nr)' * [qh; sparse(1,1,0,nn,1)];
rhs = @(t, c) rhsq2d_pod_tensor(t, c, Ar, Arp, qr, Sr);
jrhs = @(t, c) jrhsq2d_pod_tensor(t, c, Ar, Sr);

TOL = 10^(-12);
opes = odeset('AbsTol', TOL/1000, 'RelTol', TOL, 'Jacobian', jrhs, 'Mass', Mr, 'Stats', 'on');
[~, a_ode15s] = ode15s(rhs, tiempos, a0, opes);
Y_ode15s = Phi(:,1:nr)*a_ode15s' + p;

[a_bdf1,k_bdf1] = bdf1(a0, rhs, jrhs, Mr, tiempos, max_iter); % initial condition is projection on H1 -> Implicit Euler 
Y_bdf1   = Phi(:,1:nr)*a_bdf1 + p ;

[a_bdf2,k_bdf2] = bdf2(a0, rhs, jrhs, Mr, tiempos, max_iter); % initial conditions from BDF1 with delta t << bdf2 delta t
Y_bdf2   = Phi(:,1:nr)*a_bdf2 + p;

[a_bdf3,k_bdf3] = bdf3(a0, rhs, jrhs, Mr, tiempos, max_iter); % initial conditions from BDF2 with delta t << bdf3 delta t 
Y_bdf3   = Phi(:,1:nr)*a_bdf3 + p;

[a_bdf4,k_bdf4] = bdf4(a0, rhs, jrhs, Mr, tiempos, max_iter); % initial conditions from BDF3 with delta t << bdf4 delta t
Y_bdf4   = Phi(:,1:nr)*a_bdf4 + p;

[a_bdf5,k_bdf5] = bdf5(a0, rhs, jrhs, Mr, tiempos, max_iter);    % initial conditions from BDF4 with delta t << bdf5 delta t
Y_bdf5   = Phi(:,1:nr)*a_bdf5 + p;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ode15s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_ode = zeros(nn,length(tiempos)); V_ode = U_ode; 
U_ode(1:nn,:) = Y_ode15s(1:nn,:); V_ode(1:nn,:) = Y_ode15s(nn+1:end,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bdf1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_bdf1 = U_ode; V_bdf1 = U_ode; 
U_bdf1(1:nn,:) = Y_bdf1(1:nn,:); V_bdf1(1:nn,:) = Y_bdf1(nn+1:end,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bdf2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_bdf2 = U_ode; V_bdf2 = U_ode;
U_bdf2(1:nn,:) = Y_bdf2(1:nn,:); V_bdf2(1:nn,:) = Y_bdf2(nn+1:end,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bdf3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_bdf3 = U_ode; V_bdf3 = U_ode;
U_bdf3(1:nn,:) = Y_bdf3(1:nn,:); V_bdf3(1:nn,:) = Y_bdf3(nn+1:end,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bdf4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_bdf4 = U_ode; V_bdf4 = U_ode;
U_bdf4(1:nn,:) = Y_bdf4(1:nn,:); V_bdf4(1:nn,:) = Y_bdf4(nn+1:end,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bdf5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_bdf5 = U_ode; V_bdf5 = U_ode;
U_bdf5(1:nn,:) = Y_bdf5(1:nn,:); V_bdf5(1:nn,:) = Y_bdf5(nn+1:end,:);


% Compute error between bdf schemes and ode15s
[~, ~, ~, ~, errL2CombinedBdf1, errH1CombinedBdf1] = compute_errors(U_bdf1, V_bdf1, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf2, errH1CombinedBdf2] = compute_errors(U_bdf2, V_bdf2, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf3, errH1CombinedBdf3] = compute_errors(U_bdf3, V_bdf3, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf4, errH1CombinedBdf4] = compute_errors(U_bdf4, V_bdf4, U_ode, V_ode, Mh, Sh);
[~, ~, ~, ~, errL2CombinedBdf5, errH1CombinedBdf5] = compute_errors(U_bdf5, V_bdf5, U_ode, V_ode, Mh, Sh);

maxerrL2 = [max(abs(errL2CombinedBdf1)),max(abs(errL2CombinedBdf2(3:end))),max(abs(errL2CombinedBdf3(4:end))),max(abs(errL2CombinedBdf4(5:end))),max(abs(errL2CombinedBdf5(6:end)))];
maxerrH1 = [max(abs(errH1CombinedBdf1)),max(abs(errH1CombinedBdf2(3:end))),max(abs(errH1CombinedBdf3(4:end))),max(abs(errH1CombinedBdf4(5:end))),max(abs(errH1CombinedBdf5(6:end)))];

nm_iter = [k_bdf1, k_bdf2, k_bdf3, k_bdf4, k_bdf5];


save(full_name_results, 'param', 'M', 's', 'nr', 'Sh', 'Mh','Phi','Yh', 'p',...
    'Y_ode15s', 'Y_bdf1', 'Y_bdf2','Y_bdf3', 'Y_bdf4', 'Y_bdf5', ...
    'errL2CombinedBdf1','errH1CombinedBdf1',...
    'errL2CombinedBdf2','errH1CombinedBdf2',...
    'errL2CombinedBdf3','errH1CombinedBdf3',...
    'errL2CombinedBdf4','errH1CombinedBdf4',...
    'errL2CombinedBdf5','errH1CombinedBdf5', ...
    'maxerrL2', 'maxerrH1', 'nm_iter');

disp(['Results saved in file ', full_name_results,'n\.']);
