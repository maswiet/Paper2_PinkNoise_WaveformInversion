load ../data/forge2022_events.mat; load ../data/forge_picks_sensors.mat
fs=ev.fs; t=(0:size(ev.W,1)-1)'/fs-ev.pre;
bands=[10 20;20 40;40 80;80 160;160 320;320 640;640 1280];
for g=1:2
 k=find(good(:,g)); ch=(1:3)+3*(g-1); snr=zeros(numel(k),size(bands,1));
 for ib=1:size(bands,1)
  [b,a]=butter(3,bands(ib,:)/(fs/2));
  for q=1:numel(k)
   X=filtfilt(b,a,double(ev.W(:,ch,k(q)))); E=sum(X.^2,2);
   i0=t<tp(k(q),g)-0.01 & t>-0.09; i1=t>tp(k(q),g) & t<tp(k(q),g)+0.2;
   snr(q,ib)=sqrt(max(E(i1))/mean(E(i0)));
  end
 end
 fprintf('Sensor %c median amplitude SNR per band:\n',char('A'+g-1)); disp([bands(:,1)'; median(snr); mean(snr>5)])
end
