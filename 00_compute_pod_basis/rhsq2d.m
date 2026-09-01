function rhsint = rhsq2d(t, y, Tri, z, I_dir, Ahuu, Ahvv, Ahvu, fh, J, r, param)
    % Parámetros.
    nu = param(1); alpha = param(2); beta = param(3);
    
    %nt = size(Tri,1); % nt = número de triángulos
    ne = size(Tri,2); % ne = número de nodos por elemento -> si ne = 3 (lineales); ne = 6 (cuadráticos)
    nn = length(z); % nn = número de nodos para cada variable 

    % Buscamos los índices con las condiciones de contorno Dirichlet :
    %ndir = length(I_dir);
    I_all = (1:nn)';               % Todos los nodos
    I_int = setdiff(I_all, I_dir); % Nodos interiores = nodos - Dirichlet
    nint = length(I_int);

    % Antes de calcular rhs : 
    % Separamos el vector en las dos variables del problema : u,v
    u = y(1:nint); v = y(nint+1:end);
    ua = zeros(nn,1); va = ua;
    % Añadimos las condiciones Dirichlet
    ua(I_int) = u; va(I_int) = v;
    ua(I_dir) = alpha; va(I_dir) = beta/alpha;

    % % Coordenadas de los vértices de todos los triángulos.
    % x1 = z(Tri(:,1),1); x2 = z(Tri(:,2),1); x3 = z(Tri(:,3),1); 
    % y1 = z(Tri(:,1),2); y2 = z(Tri(:,2),2); y3 = z(Tri(:,3),2); 
    % % Jacobiano del cambio de variable : [xi_x xi_y; eta_x, eta_y]./J
    % J = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1); %determinante del jacobiano
    
    
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
      
    % Cálculo rhs 

    % Calculamos el término 
    % fh = [int alpha * phi1; int alpha phi2, ... , int alpha phinn] = 
    %fh = alpha*Mh*ones(nn,1);
    % Calculamos el término lineal : (con las condiciones Dirichlet
    % incluidas)
    %  fluint = (-nu*Sh(I_int,I_int) + (beta+1)*Mh(I_int,I_int))*ua(I_int) + ...
    %     (-nu*Sh(I_int,I_dir) + (beta+1)*Mh(I_int,I_dir))*gu + ...
    %     fh(I_int);
    % 
    % flvint = (-nu*Sh(I_int,I_int)*va(I_int) - Sh(I_int,I_dir)*gv) + ...
    %     beta*(Mh(I_int,I_int)*ua(I_int) + Mh(I_int,I_dir)*gu) ;
    % 
    fluint = Ahuu(I_int,I_int)*ua(I_int) + Ahuu(I_int,I_dir)*ua(I_dir) + fh(I_int);
   
    flvint = Ahvv(I_int,I_int)*va(I_int) + Ahvv(I_int,I_dir)*va(I_dir) + ...
        Ahvu(I_int,I_int)*ua(I_int) + Ahvu(I_int,I_dir)*ua(I_dir) ;
    
    % Calculamos el término no lineal nlh 
    U = ua(Tri); V = va(Tri);
    Uc = U*N; Vc = V*N;
    Nlg = (Uc.^2).*Vc; 
    NlgJ = kron(abs(J),ones(1,length(W))).*Nlg;
    Cnlg = N*diag(W)*NlgJ';
    nlh = full(sparse(Tri(:,1:ne)',ones(size(Cnlg)),Cnlg,nn,1));
 
    % Expresión final del lado derecho de la ODE en tiempo
     rhsint = [fluint + nlh(I_int) ; flvint - nlh(I_int)];
end