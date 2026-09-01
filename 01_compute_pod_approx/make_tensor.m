function [S] = make_tensor(PhiR,Tri,z,I_dir,J,l,p)
    nr = size(PhiR,2); % number of mode
    nn = length(z); % number of nodes
    
    I_all = (1:nn)';               % Todos los nodos
    I_int = setdiff(I_all, I_dir); % Nodos interiores = nodos - Dirichlet
    Iuv_int = [I_int; nn+I_int];

    Phi = zeros(size(PhiR,1),1+size(PhiR,2));
    Phi(:,1) = p; Phi(Iuv_int,2:end) = PhiR(Iuv_int,:);
    
    Phiu = Phi(1:nn,:); Phiv = Phi(nn+1:2*nn,:);
    S = zeros(nr, nr+1, nr+1, nr+1);
    for k = 1:nr+1
        for j=1:nr+1
            for i=1:k
                f = uvw(Phiu(:,i),Phiv(:,j),Phiu(:,k),Tri, z, J,l);
                % eliminamos las condiciones Dirichlet del problema a
                coef = PhiR'*[f;-f];
                S(:,i,j,k) = coef;
                S(:,k,j,i) = coef;
            end
        end
    end
end
