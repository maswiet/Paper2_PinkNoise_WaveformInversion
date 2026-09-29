load ../data/forge2022_events.mat; load ../data/tt_inversion.mat; G=struct(); q=TT(1).q; sens=[q(1:3); q(1:2) q(4)];
fs=ev.fs; t=(0:size(ev.W,1)-1)'/fs-ev.pre; [b,a]=butter(3,[40 200]/(fs/2));
for g=1:2
  k=obs(obs(:,2)==g,1); tpp=obs(obs(:,2)==g,3); ch=(1:3)+3*(g-1);
  F=zeros(numel(k),3); inc=zeros(numel(k),1); baz=inc; pol=zeros(numel(k),3);
  for i=1:numel(k)
    X=filtfilt(b,a,double(ev.W(:,ch,k(i)))); w=t>=tpp(i)-0.002 & t<tpp(i)+0.012;
    C=X(w,:)'*X(w,:); [V,D]=eig(C); [~,j]=max(diag(D)); pol(i,:)=V(:,j)';
    F(i,:)=diag(C)'/trace(C);
    d=ev.xyz(k(i),:)-sens(g,:); inc(i)=acosd(abs(d(3))/norm(d)); baz(i)=atan2d(d(1),d(2));
  end
  fprintf('sensor %c: mean P energy fraction per channel: %s ; mean incidence from vertical %.0f deg\n',char('A'+g-1),mat2str(mean(F),2),mean(inc));
  % which channel correlates with cos(incidence)?
  for c=1:3, r(c)=corr(abs(pol(:,c)),cosd(inc)); end, disp(['corr |pol_c| vs cos(inc): ' mat2str(r,2)])
  % horizontal azimuth of P in the 2 non-vertical channels vs true back-azimuth
  [~,vz]=max(r); hz=setdiff(1:3,vz); az=atan2d(pol(:,hz(1)),pol(:,hz(2)));
  sgn=sign(pol(:,vz)); az=atan2d(pol(:,hz(1)).*sgn,pol(:,hz(2)).*sgn);
  dd=wrapTo180(az-baz); fprintf('  vertical=ch%d, horizontal az offset: circ-mean %.0f deg, spread(MAD) %.0f deg\n',ch(vz),rad2deg(angle(mean(exp(1i*deg2rad(dd))))),median(abs(wrapTo180(dd-rad2deg(angle(mean(exp(1i*deg2rad(dd)))))))));
end
