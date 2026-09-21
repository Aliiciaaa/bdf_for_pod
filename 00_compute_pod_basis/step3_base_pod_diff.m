% Code for computing the POD basis
clear all; close all;

Ms = [32,64,128,256,512,1024]; 
beta = 3;

for i = 1 : length(Ms)
    tic;
    M = Ms(i);

    nombre = ['snapshots_diff_M',num2str(M),'.mat'];
    if ~isfile(nombre)
       error(['File not found: ', nombre]);
    end

    load(nombre);

    N = size(UN,2);
    nn = length(z);
    
    UNint = UN(Iuv,:);

    S2 = kron(speye(2),Sh(Iu,Iu));
    disp('Computing Cholesky factorization ...');
    [Rh,iflag,P] = chol(S2); % S2 = P*Rh'*Rh*P'
    disp('... done.');
    disp('Computing the time derivatives ...')
    % correlation matrix for POD corresponds to: (A'*A)/N
    A = Rh*(P'*UNint);
    disp('... done.');
    disp('Computing SVD ...')
    [Wl,S,Vr] = svd(A,"econ");
    disp('... done.')
    s = diag(S)/sqrt(N); % in s the singular values
    % lambda_k = s_k^2 
    
    figure();
    semilogy(s.^2, 'b', 'LineWidth', 1.5);
    title(['Correlation Matrix Eigenvalues $M=$', num2str(M)],'Interpreter','latex');
    
    PhiDiff = zeros(size(Wl));
    PhiDiff(Iuv,:) = P*(Rh\Wl); 
    time_base_pod = toc;
    save(['base_pod_diff_M',num2str(M),'.mat'], 'PhiDiff','s','UN','Zh', 'Yh','Zh','z','tp',...
        'param', 'tiempos', 'Tri', 'z', 'I_dir', 'Iu','Iv','Iuv','Ahuu', 'Ahvv', 'Ahvu', 'qh', 'Mh', 'Sh','J', 'l','p');
end
