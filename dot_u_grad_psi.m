function out = dot_u_grad_psi(u, psi, ops)
%DOT_U_GRAD_PSI Compute u · grad(psi).

[psix, psiy, psiz] = spectral_utils('grad', psi, ops);
out = u.x .* psix + u.y .* psiy + u.z .* psiz;
end
