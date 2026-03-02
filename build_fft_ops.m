function ops = build_fft_ops(n_s, L)
%BUILD_FFT_OPS Precompute wave numbers and spectral factors for periodic FFT.

k1 = (2*pi/L) * [0:(n_s/2-1), -n_s/2:-1];
[KX, KY, KZ] = ndgrid(k1, k1, k1);
K2 = KX.^2 + KY.^2 + KZ.^2;

ops.n_s = n_s;
ops.L = L;
ops.kx = KX;
ops.ky = KY;
ops.kz = KZ;
ops.k2 = K2;
end
