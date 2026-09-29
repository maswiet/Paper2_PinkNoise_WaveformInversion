load ../data/forge2022_events.mat
fs=ev.fs; t=(0:size(ev.W,1)-1)'/fs-ev.pre; [b,a]=butter(4,[5 300]/(fs/2));
figure('Position',[50 50 1400 900]);
for k=1:6
  subplot(3,2,k); hold on
  for c=1:6, x=filtfilt(b,a,double(ev.W(:,c,k))); plot(t, x/max(abs(x))*0.45 + c,'LineWidth',0.6); end
  xlim([-0.05 0.45]); title(sprintf('ev %d M=%.2f R_{cat}=%.0f m depth=%.0f m',k,ev.M(k),ev.R(k),ev.xyz(k,3)));
end
exportgraphics(gcf,'../figs/qc_forge_events.png','Resolution',100);
