function jrhs = jrhsq2d_pod_tensor(t, cp, Ar, S)
% differential of rhsq2s_pod_tensor wrt to c 
% Tensor S structure: S(m, i, j, k) where m: test mode, i,j,k: trial modes.
% Represents trilinear term: u * v * u (symmetry in i and k).

%   Writen by Alicia Garcia-Mascaraque Herrera (last modified: March 2026).
%
%   This code comes with no guarantee or warranty of any kind.
%
%   If you use this code, please cite
%
%   B. Garcia-Archilla, A. Garcia-Mascaraque and J. Novo,
%     Using BDF schemes in the temporal integration of POD-ROM methods (to appear)
%     Please check volume, pages and year of publication with Journal for proper citation.

    nr = length(cp);
    c_aug = [1; cp]; 

    % Non-linear Jacobian: ∂NL(m)/∂c(n)
    % Deriving Σᵢⱼₖ S(m,i,j,k)·cᵢ·cⱼ·cₖ leads to three terms:
    % J(m,n) = Σⱼₖ S(m,n,j,k)cⱼcₖ + Σᵢₖ S(m,i,n,k)cᵢcₖ + Σᵢⱼ S(m,i,j,n)cᵢcⱼ

    % TERM 1: Derivative w.r.t. first 'u' (index i)
    temp1 = tensorprod(S, c_aug, 4, 1);
    J1 = tensorprod(temp1, c_aug, 3, 1); % (nr x nr+1)

    % TERM 2: Derivative w.r.t. 'v' (index j)
    S_perm_j = permute(S, [1, 3, 2, 4]); % Move index j to the derivative slot
    temp2 = tensorprod(S_perm_j, c_aug, 4, 1);
    J2 = tensorprod(temp2, c_aug, 3, 1);

    % TERM 3: Derivative w.r.t. second 'u' (index k)
    S_perm_k = permute(S, [1, 4, 2, 3]); % Move index k to the derivative slot
    temp3 = tensorprod(S_perm_k, c_aug, 4, 1);
    J3 = tensorprod(temp3, c_aug, 3, 1);

    % Sum non-linear terms and remove the constant columns (index 1)
    nonlin = J1 + J2 + J3;
    jrhs = Ar + nonlin(:, 2:end);

end