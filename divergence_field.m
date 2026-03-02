function div_u = divergence_field(u, ops)
%DIVERGENCE_FIELD Compute divergence of vector field u with spectral derivatives.

[ux_x, ~, ~] = spectral_utils('grad', u.x, ops);
[~, uy_y, ~] = spectral_utils('grad', u.y, ops);
[~, ~, uz_z] = spectral_utils('grad', u.z, ops);

div_u = real(ux_x + uy_y + uz_z);
end
