function [psi1_n, psi2_n] = normalize_pair(psi1, psi2)
%NORMALIZE_PAIR Enforce |psi1|^2 + |psi2|^2 = 1 pointwise.

denom = sqrt(abs(psi1).^2 + abs(psi2).^2);
denom(denom == 0) = 1;

psi1_n = psi1 ./ denom;
psi2_n = psi2 ./ denom;
end
