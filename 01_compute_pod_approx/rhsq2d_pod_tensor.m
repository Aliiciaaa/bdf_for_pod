function rhs = rhsq2d_pod_tensor(t, c, Ar, Arp, fr, S)
% rhs for pod system using a precomputed  tensor.
% S is the tensor representing the u*v*u nonlinear term.

%   Writen by Alicia Garcia-Mascaraque Herrera (last modified: March 2026).
%
%   This code comes with no guarantee or warranty of any kind.
%
%   If you use this code, please cite
%
%   B. Garcia-Archilla, A. Garcia-Mascaraque and J. Novo,
%     Using BDF schemes in the temporal integration of POD-ROM methods (to appear)
%     Please check volume, pages and year of publication with Journal for proper citation.
  
    c_aug = [1; c];  
    
    % Compute nonlinear term: NL(m) = Σ S(m,i,j,k) * c(i) * c(j) * c(k)
    % Tensor contractions over indices k, j, and i
    nlhr = tensorprod(tensorprod(tensorprod(S, c_aug, 4, 1), c_aug, 3, 1), c_aug, 2, 1);
    
    rhs = Ar*c + Arp + nlhr + fr;
end
