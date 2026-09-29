load ../data/forge2022_events.mat
fs=ev.fs; t=(0:size(ev.W,1)-1)'/fs-ev.pre; [b,a]=butter(4,[30 1500]/(fs/2));
ne=size(ev.W,3); snr=zeros(ne,6);
X=zeros(size(ev.W),'single');
for k=1:ne, for c=1:6, x=filtfilt(b,a,double(ev.W(:,c,k))); X(:,c,k)=x;
  snr(k,c)=max(abs(x(t>0 & t<0.4)))/ (std(x(t<-0.02))+eps); end, end
disp('median/90pct SNR per channel:'); disp(median(snr)); disp(prctile(snr,90));
s=max(snr(:,[2 3 5 6]),[],2); [~,o]=sort(s,'descend'); disp(s(o(1:12))');
figure('Position',[50 50 1500 950]);
for q=1:8
  k=o(q); subplot(4,2,q); hold on
  for c=1:6, x=X(:,c,k); plot(t, x/max(abs(x))*0.45 + c,'LineWidth',0.5); end
  xlim([-0.05 0.45]); title(sprintf('ev %d M=%.2f R_{cat}=%.0f SNR=%.0f',k,ev.M(k),ev.R(k),s(k)));
end
exportgraphics(gcf,'../figs/qc_forge_events_hisnr.png','Resolution',100);
figure; k=o(1); spectrogram(double(ev.W(:,2,k)),128,120,256,fs,'yaxis'); exportgraphics(gcf,'../figs/qc_spec.png','Resolution',90);
save ../data/qc_snr.mat snr
