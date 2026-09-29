% S05  Travel-time (P onset) inversion for the deterministic model classes.
% Unknowns shared by all classes: sensor A/B positions (6) + clock offsets (2).
% Class parameters:  H: vp | G: vp, dvp/dz | VTI: vp0, eps, del (weak-anisotropy
% phase-velocity approx, straight rays). L1 misfit; AIC/BIC for comparison.
load ../data/forge2022_events.mat; load ../data/forge_picks_sensors.mat
obs = [];
for g = 1:2
    k = find(good(:,g));
    obs = [obs; k g*ones(numel(k),1) tp(k,g)]; %#ok<AGROW>
end
xe = ev.xyz(obs(:,1),:); sg = obs(:,2); to = obs(:,3); n = numel(to);
% sensors share one (vertical) well: common E,N ; depths zA,zB
p0s = [mean([S(1).p(1:2);S(2).p(1:2)]) S(1).p(3) S(2).p(3) S(1).p(4) S(2).p(4)];
% outlier screening with homogeneous model
q = fminsearch(@(q) sum(abs(to - ttmodel(q,'H',xe,sg))), [p0s 5850], optimset('MaxFunEvals',2e4,'MaxIter',2e4));
r = to - ttmodel(q,'H',xe,sg); keep = abs(r) < 5*1.4826*median(abs(r));
fprintf('outlier screening: kept %d of %d picks\n', sum(keep), numel(keep));
xe = xe(keep,:); sg = sg(keep); to = to(keep); obs = obs(keep,:); n = numel(to);

cls = {'H','G','VTI'}; x0 = {[p0s 5850], [p0s 5850 0], [p0s 5850 0 0]};
opts = optimset('MaxFunEvals',6e4,'MaxIter',6e4,'TolX',1e-4,'TolFun',1e-6);
TT = struct();
for i = 1:3
    f = @(q) sum(abs(to - ttmodel(q, cls{i}, xe, sg)));
    q = x0{i};
    for rep = 1:4, q = fminsearch(f, q, opts); end   % restarts
    r = to - ttmodel(q, cls{i}, xe, sg);
    np = numel(q);
    % Laplace likelihood: L1 -> b = mean|r|
    b = mean(abs(r)); logL = -n*log(2*b) - n;
    TT(i).cls = cls{i}; TT(i).q = q; TT(i).res = r; TT(i).mad = median(abs(r));
    TT(i).aic = 2*np - 2*logL; TT(i).bic = np*log(n) - 2*logL;
    fprintf('%-4s  MAD=%.3f ms  mean|r|=%.3f ms  AIC=%.1f  BIC=%.1f  params:', cls{i}, 1e3*TT(i).mad, 1e3*b, TT(i).aic, TT(i).bic);
    qq = q(7:end); if numel(qq)==3, qq(2:3) = 0.2*tanh(qq(2:3)); end
    fprintf(' %.4g', qq); fprintf('\n'); TT(i).par = qq;
    fprintf('      sensor E=%.0f N=%.0f zA=%.0f zB=%.0f  offsets %.2f %.2f ms\n', q(1:4), 1e3*q(5:6));
end
% residual structure: does the residual correlate between neighbouring events?
r = TT(1).res; G = sg;
for g = 1:2
    ii = find(G==g); X = xe(ii,:); rr = r(ii);
    D = squareform(pdist(X)); C = rr - mean(rr);
    edges = [0 10 20 40 80 160 320]; rho = NaN(1,numel(edges)-1);
    for e = 1:numel(edges)-1
        [a,b2] = find(triu(D>edges(e) & D<=edges(e+1),1));
        if numel(a) > 30, rho(e) = corr(C(a), C(b2)); end
    end
    fprintf('sensor %c: residual correlation vs inter-event distance bins %s m:\n   %s\n', ...
        char('A'+g-1), mat2str(edges), mat2str(rho,2));
    TT(1).rho{g} = rho; TT(1).edges = edges;
end
save ../data/tt_inversion.mat TT obs

function tt = ttmodel(q, cls, xe, sg)
        n = size(xe,1); xs = [q(1:3); q(1:2) q(4)]; off = q(5:6); q = q(7:end);
        if numel(q) >= 3, q(2:3) = 0.2*tanh(q(2:3)); end   % |eps|,|del| < 0.2
        d = xe - xs(sg,:); R = sqrt(sum(d.^2,2));
        c = abs(d(:,3))./R; s2 = 1 - c.^2;           % sin^2(theta from vertical)
        switch cls
            case 'H',   v = q(1)*ones(n,1);
            case 'G',   zm = (xe(:,3) + xs(sg,3))/2;   % mean-depth velocity (weak gradient)
                        v = q(1) + q(2)*(zm - 2400);
            case 'VTI', v = q(1)*(1 + q(3)*s2.*c.^2 + q(2)*s2.^2);
        end
        tt = off(sg)' + R./v;
        tt = tt(:);
    end
