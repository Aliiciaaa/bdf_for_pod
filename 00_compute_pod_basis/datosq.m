function [T2,z2,G2,F2]=datosq(T,z,G,F)
% Devuelve las matrices para elementos cuadr'a-
% ticos de una triangulaci'on dada por T y z
% y una lista de lados dada por G. Opcionalmente
% F es una matriz similar a G
% INPUT
%   T  matriz nt x n1  (con n1>=3) cuyas
%       componentes T(i,j), j=1,2,3 son los v'er-
%       tices del tri'angulo i-'esimo.
%   z  matriz nn x 2 con las coordenadas xy del
%       nodo j-'esimo en la fila j.
%   G matrix ng0 x n2 (con n2>=2) cuyas com-
%       ponentes G(j,1), G(j,1) son los extremos
%       del j-'esimo lado la lista dada.
%   F matriz similar a G.
% OUTPUT
%    Matrices similares a las anteriores pero correspondientes
%    a elementos cuadr'aticos.

% dimensiones de argumentos de entrada.
[nt,ntr]=size(T);  nn=size(z,1); % nt n'umero de v'ertices.
[ns,nsr]=size(G);

% Matriz de adjacencias en A, como matriz dispersa
ia=[T(:,1); T(:,2); T(:,3)]; ja=[T(:,2); T(:,3); T(:,1)];
a=ones(3*nt,1); A=sparse(ia,ja,a,nn,nn);
[Ih,Jh,ah]=find(tril(A+A')); % elementos no nulos de la parte
                                          % triangular inferior de la suma
                                          % de A y su traspuesta.
% n'umero exacto de lados es la longitud de Ih
nnew=length(Ih); 
% los nuevos nodos puntos medios de los lados.
znew=(z(Ih,:)+z(Jh,:))/2;
z2=[z;znew];
B=sparse(Ih,Jh,nn+[1:nnew]',nn,nn);
Adj=(B+B').*A; % matriz de adyacencia con el n'umero asignada
                         % a cada lado en vez de los unos que hab'ia en A.
%'Indices de puntos medios de los lados primero, segundo y
% tercero de cada tri'angulo en I1, I2, I3
I1=full(Adj((T(:,2)-1)*nn+T(:,1)));
I2=full(Adj((T(:,3)-1)*nn+T(:,2)));
I3=full(Adj((T(:,1)-1)*nn+T(:,3)));


T2=zeros(nt,ntr+3);
T2=[T(:,1:3), I1, I2, I3, T(:,4:ntr) ];
G2=zeros(ns,nsr+1);
Im=full(Adj((G(:,2)-1)*nn+G(:,1)));
G2=[G(:,1), Im, G(:,2:nsr)];

F2=[];
if nargin>3
    [ns,nsr]=size(F);
    F2=zeros(ns,nsr+1);
    Im=full(Adj((F(:,2)-1)*nn+F(:,1)));
    F2=[F(:,1), Im, F(:,2:nsr)];
end
