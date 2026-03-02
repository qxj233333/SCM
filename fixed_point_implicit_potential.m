function psi = fixed_point_implicit_potential(rhs, p_coeff, a, b, coeff_p, psi_init, params, ops)
%FIXED_POINT_IMPLICIT_POTENTIAL
% Solve (a - b*Delta)psi = rhs - coeff_p * p_coeff .* psi
% with fixed-point iterations.

psi = psi_init;
for it = 1:params.fp_max_iter
    rhs_eff = rhs - coeff_p .* p_coeff .* psi;
    psi_new = solve_helmholtz_fft(rhs_eff, a, b, ops);

    err = norm(psi_new(:) - psi(:), 2) / max(1, norm(psi_new(:), 2));
    psi = psi_new;
    if err < params.fp_tol
        break;
    end
end
end
