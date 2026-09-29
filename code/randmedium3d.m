function t = randmedium3d(sz, dx, acf, a, seed)
% RANDMEDIUM3D  Conventional random media with a correlation length, as
% competitors of the pink-noise (scale-free) field.  Same white-noise seed
% convention as PINKNOISE3D, output zero mean / unit std.
%   acf : 'exp'   exponential ACF (von Karman kappa=0.5), PSD ~ (1+k^2 a^2)^-2
%         'gauss' Gaussian ACF,                         PSD ~ exp(-k^2 a^2/4)
%   a   : correlation length (m); dx grid spacing (m)
rng(seed);
s = fftn(randn(sz));
kk = @(n) 2*pi*[0:ceil(n/2)-1, -floor(n/2):-1]/(n*dx);
[kx,ky,kz] = ndgrid(kk(sz(1)), kk(sz(2)), kk(sz(3)));
k2 = kx.^2 + ky.^2 + kz.^2;
switch acf
    case 'exp',   A = 1 ./ (1 + k2*a^2);          % sqrt of PSD
    case 'gauss', A = exp(-k2*a^2/8);
    otherwise, error('acf');
end
t = real(ifftn(s .* A));
t = (t - mean(t(:))) / std(t(:));
end
