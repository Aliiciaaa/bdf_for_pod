function nlh = uvw(u, v, ww, Tri, z, J,l)
    ne = size(Tri,2); % ne = number of nodes per element -> if ne = 3 (linear); ne = 6 (quadratic)
    nn = length(z);   % nn = number of nodes per variable

    [C,Wq] = cuad_T0; % quadrature nodes and weights for reference triangle T0
    xi = C(1,:); eta = C(2,:); 
    if l == 1
       ne = 3;
       S1 = 0.5*[1 -1 0; -1 1 0; 0 0 0];
       S12 = 0.5*[1 0 -1;-1 0 1; 0 0 0];
       S2 = S12 + S12';
       S3 = 0.5*[1 0 -1; 0 0 0; -1 0 1];
       S4 = (1/24)*[2 1 1; 1 2 1; 1 1 2];
       N = [1-xi-eta; xi; eta];
    elseif l == 2
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
        display(['The method is nor implemented for r=',num2str(r)]);
    end
      
    % non-linear term u*v*w
    U = u(Tri); V = v(Tri); W = ww(Tri);
    Uc = U*N; Vc = V*N; Wc = W*N;
    Nlg = Uc.*Vc.*Wc; 
    NlgJ = kron(abs(J),ones(1,length(Wq))).*Nlg;
    Cnlg = N*diag(Wq)*NlgJ';
    nlh = full(sparse(Tri(:,1:ne)',ones(size(Cnlg)),Cnlg,nn,1));

end