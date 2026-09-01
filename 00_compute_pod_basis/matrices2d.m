function [Mh, Sh] = matrices2d(T,z,r)

nt = size(T,1); % nt = número de triángulos
ne = size(T,2); % ne = número de nodos por elemento -> si ne = 3 (lineales); ne = 6 (cuadráticos)
nn = length(z); % nn = número de nodos

% Calculamos la matriz Ah = Ahxx + Ahyy + Mh.
% Comenzamos construyendo las matrices dispersas Ahxx, Ahyy, Mh.
Ahxx = sparse(1,1,0,nn,nn,nt*ne*ne); % definimos Ahxx como dispersa de ceros,
                                   % nn x nn, y reservmos nt*ne*ne espacio
                                   % en memoria para ella.
Ahyy = Ahxx; Mh = Ahxx;
% Completamos los huecos de la matriz.
% Primero calculamos las entradas para cada triángulo.
% si r = 1 : necesitamos los vértices de los triángulos.
% si r = 2 : necesitamos los vértices de los triángulos y los puntos
% medios.
% Coordenadas de los vértices de todos los triángulos.
x1 = z(T(:,1),1); x2 = z(T(:,2),1); x3 = z(T(:,3),1); 
y1 = z(T(:,1),2); y2 = z(T(:,2),2); y3 = z(T(:,3),2); 
% Jacobiano del cambio de variable : [xi_x xi_y; eta_x, eta_y]./J
J = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1); %determinante del jacobiano
xix = (y3 - y1)./abs(J); xiy = -(x3 - x1)./abs(J); 
etx = -(y2 - y1)./abs(J); ety = (x2 - x1)./abs(J);

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
   
ahxx = kron(J,ones(ne)).*( kron(xix.^2,S1) + kron(xix.*etx,S2) + kron(etx.^2,S3) );
ahyy = kron(J,ones(ne)).*( kron(xiy.^2,S1) + kron(xiy.*ety,S2) + kron(ety.^2,S3) );
mh = kron(J,S4);

JJ = kron(T,ones(ne,1)); Tp = T'; Tp = Tp(:); II = kron(Tp,ones(1,ne));

Ahxx = sparse(II,JJ,ahxx,nn,nn); Ahyy=sparse(II,JJ,ahyy,nn,nn);
Mh = sparse(II,JJ,mh,nn,nn);

Sh=Ahxx+Ahyy;
end