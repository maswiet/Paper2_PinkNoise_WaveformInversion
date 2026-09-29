load ../data/forge2022_events.mat; load ../data/forge_picks_sensors.mat
g=good(:,1)|good(:,2); x=ev.xyz(g,:);
disp('event E,N,Z percentiles 2/50/98:'); disp(prctile(x,[2 50 98]));
disp([S(1).p(1:3); S(2).p(1:3)]);
for p=[1.0 1.2 1.35 1.5 1.65 1.8]
  f=pinknoise3d([128 128 128],p,3); sl=[];
  for q=1:20, l=squeeze(f(randi(128),randi(128),:)); P=abs(fft(l-mean(l))).^2; k=(1:40)'; c=polyfit(log(k),log(P(k+1)),1); sl(q)=-c(1); end
  fprintf('p=%.2f -> 1D slope %.2f, std(5-pt avg)/std = %.2f\n',p,mean(sl),std(movmean(f(:),5))/std(f(:)));
end
