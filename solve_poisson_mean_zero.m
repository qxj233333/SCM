function phi = solve_poisson_mean_zero(rhs, ops)
%SOLVE_POISSON_MEAN_ZERO Solve -Delta(phi) = rhs with periodic BC and mean(phi)=0.

rhs_hat = fftn(rhs);
phi_hat = zeros(size(rhs_hat));

mask = ops.k2 > 0;
phi_hat(mask) = rhs_hat(mask) ./ ops.k2(mask);
phi_hat(~mask) = 0; % enforce mean zero

phi = real(ifftn(phi_hat));
end
