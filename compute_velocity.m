function u = compute_velocity(psi1, psi2, hbar, ops)
%COMPUTE_VELOCITY Compute u = hbar*Re(grad(conj(psi1))*i*psi1 + grad(conj(psi2))*i*psi2).

[g1x, g1y, g1z] = spectral_utils('grad', conj(psi1), ops);
[g2x, g2y, g2z] = spectral_utils('grad', conj(psi2), ops);

u.x = hbar * real(g1x .* (1i*psi1) + g2x .* (1i*psi2));
u.y = hbar * real(g1y .* (1i*psi1) + g2y .* (1i*psi2));
u.z = hbar * real(g1z .* (1i*psi1) + g2z .* (1i*psi2));
end
