function v = sgt_synth(e, M)
% SGT_SYNTH  Contract recorded strain-rate Green's tensor e [nt x 6]
% (exx eyy ezz exy exz eyz; tensor shear components) with moment tensor M.
v = e(:,1)*M(1,1) + e(:,2)*M(2,2) + e(:,3)*M(3,3) + ...
    2*(e(:,4)*M(1,2) + e(:,5)*M(1,3) + e(:,6)*M(2,3));
v = double(v);
end
