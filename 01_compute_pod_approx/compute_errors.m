function [errL2U, errH1U, errL2V, errH1V, errL2Combined, errH1Combined] = compute_errors(U, V, Upod, Vpod, Mh, Sh)
    % function that computes errors in L2 and H_0^1
    EU = U - Upod; EV = V - Vpod;
    
    % L2
    errL2U = sqrt(abs(sum(EU .* (Mh * EU), 1))); % for u
    errL2V = sqrt(abs(sum(EV .* (Mh * EV), 1))); % for v
    errL2Combined = sqrt(errL2U.^2 + errL2V.^2); % combined
    
    % H_0^1
    errH1U = sqrt(abs(sum(EU .* (Sh * EU), 1))); % for u
    errH1V = sqrt(abs(sum(EV .* (Sh * EV), 1))); % for v
    errH1Combined = sqrt(errH1U.^2 + errH1V.^2); % combined 
end