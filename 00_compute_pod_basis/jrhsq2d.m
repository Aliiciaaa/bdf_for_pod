function J = jrhsq2d(t, y, Tri, z, I_dir, Ahuu, Ahvv, Ahvu, J, r, param)
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
    % x1 = z(T(:,1),1); x2 = z(T(:,2),1); x3 = z(T(:,3),1); 
    % y1 = z(T(:,1),2); y2 = z(T(:,2),2); y3 = z(T(:,3),2); 
    % % Jacobiano del cambio de variable : [xi_x xi_y; eta_x, eta_y]./J
    % J = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1); %determinante del jacobiano
    % 
    % % Calculamos la matriz Sh, Mh.
    % [Mh, Sh] = matrices2d(T,z,r);

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

    
    % Término lineal 
    % DluDu = (-nu*Sh(I_int,I_int) + (beta+1)*Mh(I_int,I_int));
    % DlvDu = beta*Mh(I_int,I_int);
    % DlvDv = -nu*Sh(I_int,I_int);
 
    % Término no lineal 
    U = ua(Tri); V = va(Tri);
    Uc = U*N; Vc = V*N;
    % Derivada de nlu(u,v) con respecto de u: 2uv
    dnluduc1 = kron((2 * Uc .* Vc),ones(ne,1)); % cada fila repetida ne veces;
    dnluduc2 = kron(abs(J),N);
    dnluduc = dnluduc1.*dnluduc2; % Nc*diag(nlck) puestas unas debajo de otras, donde
                                % nlck es el no linal en los nodos de cuadrtura
                                % del elemento k-'esimo
    dnludu = N*diag(W)*dnluduc'; % las partial_u en cada elemento puestas unas al lado de las otras.  % derivada respecto de U
    % Derivada de nlu(u,v) con respecto de v : u^2 
    dnludvc1 = kron(Uc.^2,ones(ne,1)); % cada fila repetida ne veces;
    dnludvc = dnludvc1.*dnluduc2;
    dnludv = N*diag(W)*dnludvc';  
    
    II=kron(Tri',ones(1,ne)); Tt=Tri'; JJ=kron(Tt(:)',ones(ne,1));
    Buu = sparse(II, JJ, dnludu, nn, nn);
    Buv = sparse(II, JJ, dnludv, nn, nn);

    % Justamos las matrices por bloques.
    J = [Ahuu(I_int,I_int) + Buu(I_int,I_int), Buv(I_int,I_int); Ahvu(I_int,I_int)-Buu(I_int,I_int), Ahvv(I_int,I_int)-Buv(I_int,I_int)];
end