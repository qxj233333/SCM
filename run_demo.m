% run_demo.m
% Demo script for SI-BDF1 + SI-BDF2 solver.
clear; clc;

params.L = 2*pi;
params.n_s = 16;
params.T = 0.2;
params.n_t = 20;
params.n_tE = 5;          % SI-BDF1 sub-steps used to generate psi^1
params.hbar = 1.0;
params.epsilon = 0.2;
params.fp_max_iter = 20;  % fixed-point iterations for implicit potential term
params.fp_tol = 1e-8;

% build grid
x = (0:params.n_s-1) * (params.L/params.n_s);
[X, Y, Z] = ndgrid(x, x, x);

% sample normalized initial condition
psi1_0 = exp(1i * (sin(X) + 0.2*cos(Y)));
psi2_0 = 0.3 * exp(1i * (cos(Y) + 0.1*sin(Z)));

norm0 = sqrt(abs(psi1_0).^2 + abs(psi2_0).^2);
psi1_0 = psi1_0 ./ norm0;
psi2_0 = psi2_0 ./ norm0;

[psi1_end, psi2_end, hist] = si_bdf2_solver(psi1_0, psi2_0, params);

fprintf('Done. Final time = %.4f\n', params.T);
fprintf('Final mean(|psi1|^2+|psi2|^2) = %.12f\n', mean(abs(psi1_end).^2 + abs(psi2_end).^2, 'all'));
fprintf('Stored %d SI-BDF2 steps.\n', numel(hist.p));
