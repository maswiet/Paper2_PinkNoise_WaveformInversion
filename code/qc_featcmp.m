G=forge_setup(); nm={'clevel2','cdecay2','swidth','spd','sp'};
Df=load('../data/forge_features.mat'); mods={'H','L','P1','P3'};
for g=1:2, for ib=1:2
  fprintf('\nsensor %c band %d-%d\n',char('A'+g-1),Df.bands(ib,:));
  F=Df.D(g,ib).F; fprintf('%-5s','data'); for j=1:5, v=[F.(nm{j})]; fprintf(' %s %7.3g', nm{j}, median(v,'omitnan')); end; fprintf('\n');
  for m=1:numel(mods), L=load(sprintf('../results/feat_%s_forge.mat',mods{m})); F=L.M(g,ib).F; [~,im]=intersect(L.M(g,ib).k,Df.D(g,ib).k);
    fprintf('%-5s',mods{m}); for j=1:5, v=[F(im).(nm{j})]; fprintf(' %s %7.3g', nm{j}, median(v,'omitnan')); end; fprintf('\n'); end
end, end
% envelope figure
figure('Position',[50 50 1200 700]); lag=(-0.08:1e-3:0.12)'; 
for g=1:2, for ib=1:2
 subplot(2,2,2*(g-1)+ib); E=[Df.D(g,ib).F.envS]; semilogy(lag,median(E,2,'omitnan'),'k','LineWidth',2); hold on
 for m=1:numel(mods), L=load(sprintf('../results/feat_%s_forge.mat',mods{m})); [~,im]=intersect(L.M(g,ib).k,Df.D(g,ib).k); E=[L.M(g,ib).F(im).envS]; semilogy(lag,median(E,2,'omitnan'),'LineWidth',1.3); end
 legend(['data' mods]); title(sprintf('sensor %c %d-%d Hz',char('A'+g-1),Df.bands(ib,:))); ylim([1e-3 2]); grid on; xlabel('t - t_S (s)');
end, end
exportgraphics(gcf,'../figs/qc_env_cmp.png','Resolution',110);
