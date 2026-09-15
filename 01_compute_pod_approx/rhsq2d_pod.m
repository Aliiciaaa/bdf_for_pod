function rhsint = rhsq2d_pod(t, c, PhiR, Tri, z, I_dir, J,r,Ar, Arp, fr, p, param)
    % Parámetros
    %nt = size(Tri,1); % nt = número de triángulos
    ne = size(Tri,2); % ne = número de nodos por elemento -> si ne = 3 (lineales); ne = 6 (cuadráticos)
    nn = length(z); % nn = número de nodos para cada variable 
    
    nu = param(1); alpha = param(2); beta = param(3);

    % Buscamos los índices con las condiciones de contorno Dirichlet :
    %ndir = length(I_dir);
    I_all = (1:nn)';               % Todos los nodos
    I_int = setdiff(I_all, I_dir); % Nodos interiores = nodos - Dirichlet
    nint = length(I_int);

    % Antes de 
    % calcular rhs : 
    % Separamos el vector en las dos variables del problema : u,v
    yh = PhiR*c; zh = yh + p;
    ua = zh(1:nn); va = zh(nn+1:end);
    % u = y(1:nint); v = y(nint+1:end);
    % ua = zeros(nn,1); va = ua;
    % % Añadimos las condiciones Dirichlet
    % ua(I_int) = u; va(I_int) = v;
    % ua(I_dir) = alpha; va(I_dir) = beta/alpha;

    % Definimos los nodos de cuadratura, las funciones base en triángulo
    % de referencia.
    [C,W] = cuad_T0; % nodos y pesos de cuadratura en T0.
    xi = C(1,:); eta = C(2,:); % para calcular el término f 
    if r == 1
       ne = 3;
       S1 = 0.5*[1 -1 0; -1 1 0; 0 0 0];
       S12 = 0.5*[1 0 -1;-1 0 1; 0 0 0];
       S2 = S12 + S12';
       S3 = 0.5*[1 0 -1; 0 0 0; -1 0 1];
       S4 = (1/24)*[2 1 1; 1 2 1; 1 1 2];
       N = [1-xi-eta; xi; eta];
    elseif r == 2
       ne = 6;  % ne = nodos por elemento
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
        display(['El método no está implementado para r=',num2str(r)]);
    end
      
    % Cálculo rhs POD
    % Calculamos el término no lineal nlh 
    U = ua(Tri); V = va(Tri);
    Uc = U*N; Vc = V*N;
    Nlg = (Uc.^2).*Vc; 
    NlgJ = kron(abs(J),ones(1,length(W))).*Nlg;
    Cnlg = N*diag(W)*NlgJ';
    nlh = full(sparse(Tri(:,1:ne)',ones(size(Cnlg)),Cnlg,nn,1));
    nlhr = PhiR'*[nlh;-nlh];
    
    % Expresión fiWnal del lado derecho de la ODE en tiempo
     rhsint = Ar*c + Arp + nlhr + fr;
end