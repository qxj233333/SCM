function p_eps = compute_pressure(psi1, psi2, u, params, ops)
%COMPUTE_PRESSURE Compute p_{epsilon,psi} by Poisson solve with mean-zero constraint.
%
% -Delta p_eps = div(-hbar * v_eps)

hbar = params.hbar;
eps2 = params.epsilon^2;

% convective-like part
u2 = u.x.^2 + u.y.^2 + u.z.^2;
[g_u2x, g_u2y, g_u2z] = spectral_utils('grad', u2, ops);

grad_u2.x = g_u2x;
grad_u2.y = g_u2y;
grad_u2.z = g_u2z;

% (u·grad)u
ugrad_u = advect_vector(u, u, ops);

term1.x = -(1-eps2) * ( (1/(2*hbar))*grad_u2.x + (1/hbar)*ugrad_u.x );
term1.y = -(1-eps2) * ( (1/(2*hbar))*grad_u2.y + (1/hbar)*ugrad_u.y );
term1.z = -(1-eps2) * ( (1/(2*hbar))*grad_u2.z + (1/hbar)*ugrad_u.z );

% quantum-like correction using psi1, psi2
q = quantum_vector_term(psi1, ops) + quantum_vector_term(psi2, ops);
term2.x = eps2 * (hbar/2) * real(q.x);
term2.y = eps2 * (hbar/2) * real(q.y);
term2.z = eps2 * (hbar/2) * real(q.z);

v.x = term1.x + term2.x;
v.y = term1.y + term2.y;
v.z = term1.z + term2.z;

rhs = divergence_field(struct('x', -hbar*v.x, 'y', -hbar*v.y, 'z', -hbar*v.z), ops);
p_eps = solve_poisson_mean_zero(rhs, ops);

% remove mean explicitly to avoid drift
p_eps = p_eps - mean(p_eps, 'all');
end

function out = quantum_vector_term(psi, ops)
lap_conj = spectral_utils('lap', conj(psi), ops);
[g_lapx, g_lapy, g_lapz] = spectral_utils('grad', lap_conj, ops);
[g_psix, g_psiy, g_psiz] = spectral_utils('grad', psi, ops);

out.x = g_lapx .* psi - lap_conj .* g_psix;
out.y = g_lapy .* psi - lap_conj .* g_psiy;
out.z = g_lapz .* psi - lap_conj .* g_psiz;
end

function out = advect_vector(u, v, ops)
[dvxx, dvxy, dvxz] = spectral_utils('grad', v.x, ops);
[dvyx, dvyy, dvyz] = spectral_utils('grad', v.y, ops);
[dvzx, dvzy, dvzz] = spectral_utils('grad', v.z, ops);

out.x = u.x.*dvxx + u.y.*dvxy + u.z.*dvxz;
out.y = u.x.*dvyx + u.y.*dvyy + u.z.*dvyz;
out.z = u.x.*dvzx + u.y.*dvzy + u.z.*dvzz;
end
