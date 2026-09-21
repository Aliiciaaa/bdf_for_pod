function jrhsint = jrhsq2d(t, y, Tri, z, I_dir, Ahuu, Ahvv, Ahvu, J, r, param)
    % Jacobian of right-hand side of the FEM system for interior nodes
    nu = param(1); alpha = param(2); beta = param(3);
    
    %nt = size(Tri,1); % nt = number of triangles 
    ne = size(Tri,2);  % ne = number of nodes per element -> if ne = 3 (linear); ne = 6 (quadratic)
    nn = length(z);    % nn = number of nodes per variable

    I_all = (1:nn)';               % All nodes 
    I_int = setdiff(I_all, I_dir); % Interior nodes = All nodes - Dirichlet nodes
    nint = length(I_int);
    
    u = y(1:nint); v = y(nint+1:end);
    ua = zeros(nn,1); va = ua;
    % add Dirichlet conditions for each variable
    ua(I_int) = u; va(I_int) = v;
    ua(I_dir) = alpha; va(I_dir) = beta/alpha;

    [C,W] = cuad_T0; % quadrature nodes and weights at reference triangle T0
    xi = C(1,:); eta = C(2,:); 
    if r == 1
       ne = 3;
       S1 = 0.5*[1 -1 0; -1 1 0; 0 0 0];
       S12 = 0.5*[1 0 -1;-1 0 1; 0 0 0];
       S2 = S12 + S12';
       S3 = 0.5*[1 0 -1; 0 0 0; -1 0 1];
       S4 = (1/24)*[2 1 1; 1 2 1; 1 1 2];
       N = [1-xi-eta; xi; eta];
    elseif r == 2
       ne = 6; 
       S1=[ 3 1 0 -4 0 0; 1 3 0 -4 0 0; 0 0 0 0 0 0; -4 -4 0 8 0 0; ...
            0 0 0 0 8 -8; 0 0 0 0 -8 8 ]/6;
       S2=[ 6 1 1 -4 0 -4; 1 0 -1 -4 4 0;  1 -1 0 0 4 -4; -4 -4 0 8 -8 8; ...
            0 4 4 -8 8 -8; -4 0 -4 8 -8 8 ]/6;
       S3=[ 3 0 1 0 0 -4; 0 0 0 0 0 0; 1 0 3 0 0 -4; 0 0 0 8 -8 0; ...
            0 0 0 -8 8 0; -4 0 -4 0 0 8]/6;
       S12=[3 0 1 0 0 -4; 1 0 -1 -4 4 0; 0 0 0 0 0 0; -4 0 0 4 -4 4; ...
            0 0 4 -4 4 -4; 0 0 -4 4 -4 4]/6;
       S4=[ 6 -1 -1  0 -4  0; -1  6 -1  0  0 -4; -1 -1  6 -4  0  0; ...
            0  0 -4 32 16 16; -4  0  0 16 32 16; 0 -4  0 16 16 32]/360;
       N = [ (1 - xi - eta).*(1 - 2*xi - 2*eta);
              xi.*(2*xi - 1);
              eta.*(2*eta - 1);
              4*xi.*(1 - xi - eta);
              4*eta.*xi;
              4*eta.*(1 - xi - eta)];
    else
        display(['The method is not implemented for r=',num2str(r)]);
    end
 
    % non-linear term 
    U = ua(Tri); V = va(Tri);
    Uc = U*N; Vc = V*N;
    % d nlu(u,v)/du = 2uv
    dnluduc1 = kron((2 * Uc .* Vc),ones(ne,1)); % repeat each row ne times
    dnluduc2 = kron(abs(J),N);
    dnluduc = dnluduc1.*dnluduc2; % Nc*diag(nlck) one below another and nlck is the non-linear term in quadrature nodes for the k-th element

    dnludu = N*diag(W)*dnluduc'; % d nlu/du for each element 

    % d nlu(u,v)/dv = u^2 
    dnludvc1 = kron(Uc.^2,ones(ne,1)); % repeat each row ne times
    dnludvc = dnludvc1.*dnluduc2;

    dnludv = N*diag(W)*dnludvc'; % d nlu/dv for each element 
    
    II = kron(Tri',ones(1,ne)); Tt=Tri'; JJ=kron(Tt(:)',ones(ne,1));
    Buu = sparse(II, JJ, dnludu, nn, nn);
    Buv = sparse(II, JJ, dnludv, nn, nn);
    % non-linear term corresponds to [Buu, Buv; -Buu, -Buv]

    jrhsint = [Ahuu(I_int,I_int) + Buu(I_int,I_int), Buv(I_int,I_int); Ahvu(I_int,I_int)-Buu(I_int,I_int), Ahvv(I_int,I_int)-Buv(I_int,I_int)];
end
