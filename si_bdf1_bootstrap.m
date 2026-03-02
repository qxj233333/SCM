function [psi1, psi2] = si_bdf1_bootstrap(psi1_0, psi2_0, dt_target, params, ops)
%SI_BDF1_BOOTSTRAP Use SI-BDF1 sub-steps to produce psi^1 for SI-BDF2.

dtE = dt_target / params.n_tE;
psi1 = psi1_0;
psi2 = psi2_0;

u = compute_velocity(psi1, psi2, params.hbar, ops);
p = compute_pressure(psi1, psi2, u, params, ops);

for n = 1:params.n_tE
    conv1 = dot_u_grad_psi(u, psi1, ops);
    conv2 = dot_u_grad_psi(u, psi2, ops);

    rhs1 = psi1 - (1-params.epsilon^2)*dtE*conv1;
    rhs2 = psi2 - (1-params.epsilon^2)*dtE*conv2;

    % (1 - i*eps^2*hbar/2*dtE*Delta) psi* = rhs - i*(1/hbar)*dtE*p.*psi*
    a = 1;
    b = 1i * params.epsilon^2 * params.hbar * dtE / 2;
    coeff_p = 1i * (1/params.hbar) * dtE;

    psi1_ast = fixed_point_implicit_potential(rhs1, p, a, b, coeff_p, psi1, params, ops);
    psi2_ast = fixed_point_implicit_potential(rhs2, p, a, b, coeff_p, psi2, params, ops);

    [psi1_ast2, psi2_ast2] = normalize_pair(psi1_ast, psi2_ast);

    u_ast2 = compute_velocity(psi1_ast2, psi2_ast2, params.hbar, ops);
    phi = solve_poisson_mean_zero(divergence_field(u_ast2, ops), ops);

    phase = exp(1i * phi / params.hbar);
    psi1 = phase .* psi1_ast2;
    psi2 = phase .* psi2_ast2;

    u = compute_velocity(psi1, psi2, params.hbar, ops);
    p = compute_pressure(psi1, psi2, u, params, ops);
end
end
