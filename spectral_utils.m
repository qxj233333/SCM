function varargout = spectral_utils(action, field, ops)
%SPECTRAL_UTILS Basic spectral derivative utilities (periodic domain).

switch action
    case 'grad'
        fhat = fftn(field);
        gx = ifftn(1i * ops.kx .* fhat);
        gy = ifftn(1i * ops.ky .* fhat);
        gz = ifftn(1i * ops.kz .* fhat);
        varargout = {gx, gy, gz};

    case 'lap'
        fhat = fftn(field);
        lap = ifftn(-(ops.k2) .* fhat);
        varargout = {lap};

    otherwise
        error('Unknown action in spectral_utils: %s', action);
end
end
