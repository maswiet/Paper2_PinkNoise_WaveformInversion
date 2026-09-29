function t = pinknoise3d(sz, p, seed)
% PINKNOISE3D  Spatially correlated 3D "pink-noise" field, GFI recipe
% (Web Session 4, slide 7): Gaussian white noise, Fourier amplitude
% filtered by 1/(1+|k|)^p, back-transformed.  Output normalised to
% zero mean, unit std.  p = 1.35 gives 1D well-log spectra S(k) ~ 1/k^~1.
%
%   sz   : [nx ny nz]
%   p    : 3D filter exponent (p=0 -> white noise)
%   seed : rng seed (for reproducible realisations)
if nargin > 2, rng(seed); end
r = randn(sz);
s = fftn(r);
[kx,ky,kz] = ndgrid(fk(sz(1)), fk(sz(2)), fk(sz(3)));
R = 1 + sqrt(kx.^2 + ky.^2 + kz.^2);
t = real(ifftn(s ./ R.^p));
t = (t - mean(t(:))) / std(t(:));
end

function k = fk(n)
% integer wavenumbers in FFT order, 0..n/2-1, -n/2..-1
k = [0:ceil(n/2)-1, -floor(n/2):-1];
end
