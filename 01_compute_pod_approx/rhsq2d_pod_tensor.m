function rhs = rhsq2d_pod_tensor(t, c, Ar, Arp, fr, S)
    % rhs for pod system using a precomputed  tensor.
    % S is the tensor representing the u*v*u nonlinear term.
    c_aug = [1; c];  
    
    % NL(m) = Σ S(m,i,j,k) * c(i) * c(j) * c(k)
    nlhr = tensorprod(tensorprod(tensorprod(S, c_aug, 4, 1), c_aug, 3, 1), c_aug, 2, 1);
    
    rhs = Ar*c + Arp + nlhr + fr;
end
