function x = solve_helmholtz_fft(rhs, a, b, ops)
%SOLVE_HELMHOLTZ_FFT Solve (a - b*Delta)x = rhs with periodic BC by FFT.

rhs_hat = fftn(rhs);
denom = a + b * ops.k2;

x_hat = rhs_hat ./ denom;
x = ifftn(x_hat);
end
