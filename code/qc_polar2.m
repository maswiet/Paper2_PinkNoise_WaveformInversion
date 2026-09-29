load ../data/forge2022_events.mat; load ../data/tt_inversion.mat; q=TT(1).q; sens=[q(1:3); q(1:2) q(4)];
fs=ev.fs; t=(0:size(ev.W,1)-1)'/fs-ev.pre; [b,a]=butter(3,[40 200]/(fs/2));
for g=1:2
  k=obs(obs(:,2)==g,1); tpp=obs(obs(:,2)==g,3); ch=(1:3)+3*(g-1);
  az=zeros(numel(k),1); baz=az; lin=az;
  for i=1:numel(k)
    X=filtfilt(b,a,double(ev.W(:,ch(2:3),k(i)))); w=t>=tpp(i)-0.002 & t<tpp(i)+0.010;
    C=X(w,:)'*X(w,:); [V,D]=eig(C); [dm,j]=max(diag(D)); lin(i)=1-min(diag(D))/dm;
    az(i)=atan2d(V(1,j),V(2,j));   % axial (180 deg ambiguity)
    d=ev.xyz(k(i),:)-sens(g,:); baz(i)=atan2d(d(1),d(2));
  end
  dd=mod(az-baz,180); m=rad2deg(angle(mean(exp(2i*deg2rad(dd)))))/2;
  fprintf('sensor %c (ch%d,ch%d horizontals): P rectilinearity %.2f, axial az offset %.0f deg, spread(MAD) %.0f deg, event baz range %.0f-%.0f\n', ...
    char('A'+g-1),ch(2),ch(3),median(lin),mod(m,180),median(abs(mod(dd-m+90,180)-90)),prctile(baz,5),prctile(baz,95));
end
