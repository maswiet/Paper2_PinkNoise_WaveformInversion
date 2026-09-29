function M = dc_moment(strike, dip, rake)
% DC_MOMENT  Unit double-couple moment tensor (Aki & Richards 1980, Box 4.4),
% in the x=North, y=East, z=Down convention, returned in the FD frame
% x=East, y=North, z=Down.
f = deg2rad(strike); d = deg2rad(dip); l = deg2rad(rake);
Mxx = -(sin(d)*cos(l)*sin(2*f) + sin(2*d)*sin(l)*sin(f)^2);
Mxy =  (sin(d)*cos(l)*cos(2*f) + 0.5*sin(2*d)*sin(l)*sin(2*f));
Mxz = -(cos(d)*cos(l)*cos(f)   + cos(2*d)*sin(l)*sin(f));
Myy =  (sin(d)*cos(l)*sin(2*f) - sin(2*d)*sin(l)*cos(f)^2);
Myz = -(cos(d)*cos(l)*sin(f)   - cos(2*d)*sin(l)*cos(f));
Mzz =  sin(2*d)*sin(l);
Mar = [Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz];   % (N,E,D)
P = [0 1 0; 1 0 0; 0 0 1];                        % swap N<->E
M = P*Mar*P';
end
