function [psi1, psi2, history] = si_bdf2_solver(psi1_0, psi2_0, params)
%SI_BDF2_SOLVER Semi-implicit BDF2 for the constrained (psi1, psi2) system.
%
% Inputs
%   psi1_0, psi2_0 : initial complex fields (n_s x n_s x n_s)
%   params         : struct with required fields:
%                    L, n_s, T, n_t, n_tE, hbar, epsilon,
%                    fp_max_iter, fp_tol
%
% Outputs
%   psi1, psi2     : final fields at time T
%   history        : struct containing intermediate pressure and velocity

validateattributes(psi1_0, {'numeric'}, {'nonempty'});
validateattributes(psi2_0, {'numeric'}, {'size', size(psi1_0)});

dt = params.T / params.n_t;

% Build spectral wave numbers and helper operators
ops = build_fft_ops(params.n_s, params.L);

% --- SI-BDF1 bootstrap: obtain psi^1 from psi^0 ---
[psi1_1, psi2_1] = si_bdf1_bootstrap(psi1_0, psi2_0, dt, params, ops);

psi1_prev = psi1_0;
psi2_prev = psi2_0;
psi1 = psi1_1;
psi2 = psi2_1;

u_prev = compute_velocity(psi1_prev, psi2_prev, params.hbar, ops);
u = compute_velocity(psi1, psi2, params.hbar, ops);

p_prev = compute_pressure(psi1_prev, psi2_prev, u_prev, params, ops);
p = compute_pressure(psi1, psi2, u, params, ops);

history.p = cell(params.n_t + 1, 1);
history.u = cell(params.n_t + 1, 1);
history.p{1} = p_prev;
history.p{2} = p;
history.u{1} = u_prev;
history.u{2} = u;

for n = 1:(params.n_t - 1)
    % AB2 extrapolation for convection and pressure coefficient
    conv1_ext = 2 * dot_u_grad_psi(u, psi1, ops) - dot_u_grad_psi(u_prev, psi1_prev, ops);
    conv2_ext = 2 * dot_u_grad_psi(u, psi2, ops) - dot_u_grad_psi(u_prev, psi2_prev, ops);
    p_ext = 2 * p - p_prev;

    rhs1 = 4*psi1 - psi1_prev - 2*(1-params.epsilon^2)*dt*conv1_ext;
    rhs2 = 4*psi2 - psi2_prev - 2*(1-params.epsilon^2)*dt*conv2_ext;

    % Solve predictor with fixed-point iteration:
    % (3 - i*eps^2*hbar*dt*Delta) psi* = rhs - i*(2/hbar)*dt*p_ext.*psi*
    a = 3;
    b = 1i * params.epsilon^2 * params.hbar * dt;
    coeff_p = 1i * (2/params.hbar) * dt;

    psi1_ast = fixed_point_implicit_potential(rhs1, p_ext, a, b, coeff_p, psi1, params, ops);
    psi2_ast = fixed_point_implicit_potential(rhs2, p_ext, a, b, coeff_p, psi2, params, ops);

    % Projection 1: normalization
    [psi1_ast2, psi2_ast2] = normalize_pair(psi1_ast, psi2_ast);

    % Projection 2: incompressibility via phase correction
    u_ast2 = compute_velocity(psi1_ast2, psi2_ast2, params.hbar, ops);
    phi = solve_poisson_mean_zero(divergence_field(u_ast2, ops), ops);

    phase = exp(1i * phi / params.hbar);
    psi1_next = phase .* psi1_ast2;
    psi2_next = phase .* psi2_ast2;

    % Update u and pressure
    u_next = compute_velocity(psi1_next, psi2_next, params.hbar, ops);
    p_next = compute_pressure(psi1_next, psi2_next, u_next, params, ops);

    % Shift time levels
    psi1_prev = psi1; psi2_prev = psi2;
    psi1 = psi1_next; psi2 = psi2_next;

    u_prev = u; u = u_next;
    p_prev = p; p = p_next;

    history.p{n+2} = p;
    history.u{n+2} = u;
end
end
