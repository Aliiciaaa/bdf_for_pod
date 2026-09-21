function rhsint = rhsq2d(t, y, Tri, z, I_dir, Ahuu, Ahvv, Ahvu, fh, J, r, param)
    % Right-hand side of the FEM system for interior nodes.
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

    % option 1. computing the linear part using stiffness and mass matrix
    % from FEM system for interior nodes
    % fh = [int alpha * phi1; int alpha phi2, ... , int alpha phinn] = 
    % fh = alpha*Mh*ones(nn,1);
    % 
    % fluint = (-nu*Sh(I_int,I_int) + (beta+1)*Mh(I_int,I_int))*ua(I_int) + ...
    %     (-nu*Sh(I_int,I_dir) + (beta+1)*Mh(I_int,I_dir))*gu + ...
    %     fh(I_int);
    % 
    % flvint = (-nu*Sh(I_int,I_int)*va(I_int) - Sh(I_int,I_dir)*gv) + ...
    %     beta*(Mh(I_int,I_int)*ua(I_int) + Mh(I_int,I_dir)*gu) ;
    % 
    % option 2. computing the linear part using Ah matricess (usefull for
    % the POD system)
    fluint = Ahuu(I_int,I_int)*ua(I_int) + Ahuu(I_int,I_dir)*ua(I_dir) + fh(I_int);
   
    flvint = Ahvv(I_int,I_int)*va(I_int) + Ahvv(I_int,I_dir)*va(I_dir) + ...
        Ahvu(I_int,I_int)*ua(I_int) + Ahvu(I_int,I_dir)*ua(I_dir) ;
    
    % non-linear term 
    U = ua(Tri); V = va(Tri);
    Uc = U*N; Vc = V*N;
    Nlg = (Uc.^2).*Vc; 
    NlgJ = kron(abs(J),ones(1,length(W))).*Nlg;
    Cnlg = N*diag(W)*NlgJ';
    nlh = full(sparse(Tri(:,1:ne)',ones(size(Cnlg)),Cnlg,nn,1));
 
    rhsint = [fluint + nlh(I_int) ; flvint - nlh(I_int)];
end
