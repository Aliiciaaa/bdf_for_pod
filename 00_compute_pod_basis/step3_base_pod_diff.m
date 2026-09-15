% Code for computing the POD basis
clear all; close all;

Ms = [32,64,128,256,512]; 
beta = 3;

for i = 1 : length(Ms)
    tic;
    M = Ms(i);

    nombre = ['snapshots_diff_M',num2str(M),'.mat'];
    if ~isfile(nombre)
       error(['Archivo no encontrado: ', nombre]);
    end

    load(nombre);

    N = size(UN,2);
    nn = length(z);
    
    UNint = UN(Iuv,:);
    % Producto escalar en H1 elementos finitos sobre el espacio.
    S2 = kron(speye(2),Sh(Iu,Iu));
    disp('Computing Cholesky factorization ...');
    [Rh,iflag,P] = chol(S2); % S2 = P*Rh'*Rh*P'
    disp('... done.');
    disp('computing the time derivatives ...')
    % La matriz de correlación del método POD la denotamos como (A'*A)/N
    A = Rh*(P'*UNint);
    % Calculamos la base POD
    disp('Computing SVD ...')
    [Wl,S,Vr] = svd(A,"econ");
    disp('... done.')
    s = diag(S)/sqrt(N); % in s the singular values
    % lambda_k = s_k^2 
    figure();
    semilogy(s.^2, 'b', 'LineWidth', 1.5);
    title('Autovalores de la matriz de correlación','Interpreter','latex');
    
    PhiDiff = zeros(size(Wl));
    PhiDiff(Iuv,:) = P*(Rh\Wl); 
    time_base_pod = toc;
    save(['base_pod_diff_M',num2str(M),'.mat'], 'PhiDiff','s','UN','Zh', 'Yh','Zh','z','tp',...
        'param', 'tiempos', 'Tri', 'z', 'I_dir', 'Iu','Iv','Iuv','Ahuu', 'Ahvv', 'Ahvu', 'qh', 'Mh', 'Sh','J', 'l','p');
end
