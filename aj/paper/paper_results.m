function paper_results(outDir)
%PAPER_RESULTS  Navigation-level results of the 70-run campaign (R1 Section VII).
%   paper_results(outDir)
%   Reads the per-item results of run_campaign (ajPath('campaign','results'),
%   <scen>_<profile>_s<seed>.mat) and writes to outDir (default: <repo>\paper\R1)
%     paper_results.txt  the values of Tables III and IV and of Sections VII-B/C
%                        (median and 95th-percentile rms error, time with fewer
%                        than four tracking channels, mean tracking channels,
%                        S1 J/S tolerance, paired ratios, seed-3 aiding ablation)
%     fig_campaign.pdf   Fig. 9 (rms 3-D position error per run)
%   The horizontal jitter of the markers in Fig. 9 is random (seeded here).
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo, 'aj'));
if nargin < 1, outDir = fullfile(repo, 'paper', 'R1'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
txt = fullfile(outDir, 'paper_results.txt');
out = fopen(txt, 'w');  cl = onCleanup(@() fclose(out));
P = @(varargin) fprintf(out, varargin{:});
F = dir(fullfile(ajPath('campaign', 'results'), 'S*.mat'));      % excludes noins_*.mat
cfg = {'single','pi','pi_adn','mvdr','mvdr_adn','pi_a','pi_dn','pi_ab','mvdr_ab','pi_ff'};
fld = {'rmsJam','lostFrac','chTrackMean','tFail','jsTol','tRecover','maxJam','coastMax'};
n = numel(F);
D = struct();  for c = cfg, for f = fld, D.(c{1}).(f{1}) = NaN(n,1); end, D.(c{1}).ctRms = NaN(n,1); D.(c{1}).ctJs = NaN(n,1); end
sc = cell(n,1); pr = cell(n,1); sd = NaN(n,1);
for i = 1:n
  R = load(fullfile(F(i).folder, F(i).name));
  sc{i} = R.scen; pr{i} = R.profile; sd(i) = R.seed;
  for c = cfg
    if ~isfield(R, c{1}), continue; end
    u = R.(c{1}).UT;
    for f = fld, if isfield(u, f{1}) && ~isempty(u.(f{1})), D.(c{1}).(f{1})(i) = double(u.(f{1})); end, end
    k = R.(c{1}).CT;
    if isfield(k,'rmsJam'), D.(c{1}).ctRms(i) = k.rmsJam; end
    if isfield(k,'jsTol') && ~isempty(k.jsTol), D.(c{1}).ctJs(i) = k.jsTol; end
  end
end
P('items %d\n', n);
for s = {'S1','S2','S3'}, P('%s: %d items (%s)\n', s{1}, sum(strcmp(sc,s{1})), strjoin(unique(pr(strcmp(sc,s{1})))', ',')); end
p95 = @(x) prctile(x(~isnan(x)), 95);
main = {'single','pi','pi_adn','mvdr','mvdr_adn'};
P('\n== Table III: per scenario, UT rms median / p95 [m], lost %%, track mean, failures, recover median, CT rms median\n');
for s = {'S1','S2','S3','all'}
  k = true(n,1); if ~strcmp(s{1},'all'), k = strcmp(sc,s{1}); end
  for c = main
    d = D.(c{1});
    P('%-3s %-9s n=%2d med %7.3f p95 %7.3f lost %5.1f track %4.2f fail %2d/%2d rec %5.2f CT %8.2f\n', s{1}, c{1}, sum(k), ...
      median(d.rmsJam(k),'omitnan'), p95(d.rmsJam(k)), 100*mean(d.lostFrac(k),'omitnan'), mean(d.chTrackMean(k),'omitnan'), ...
      sum(isfinite(d.tFail(k))), sum(k), median(d.tRecover(k),'omitnan'), median(d.ctRms(k),'omitnan'));
  end
end
P('\n== per scenario/profile median UT rms\n');
for s = {'S1','S2','S3'}
  for p = unique(pr(strcmp(sc,s{1})))'
    k = strcmp(sc,s{1}) & strcmp(pr,p{1});
    P('%s %-16s n=%2d', s{1}, p{1}, sum(k));
    for c = main, P(' %s %.3f', c{1}, median(D.(c{1}).rmsJam(k),'omitnan')); end
    P(' | lost pi %.1f adn %.1f mvdr %.1f\n', 100*mean(D.pi.lostFrac(k)), 100*mean(D.pi_adn.lostFrac(k)), 100*mean(D.mvdr.lostFrac(k)));
  end
end
P('\n== S1 J/S tolerance (UT; Inf = no failure up to 100 dB), CT tolerance\n');
k = strcmp(sc,'S1');
for c = main
  j = D.(c{1}).jsTol(k);  jc = D.(c{1}).ctJs(k);
  P('%-9s UT med %6.1f min %6.1f noFail %2d/%2d | CT med %6.1f noFail %2d\n', c{1}, median(j,'omitnan'), min(j), sum(isinf(j)), numel(j), median(jc,'omitnan'), sum(isinf(jc)));
end
P('\n== Table IV: paired ratios (median, [p25 p75]), runs in which the first is better\n');
pairs = {'mvdr','pi'; 'pi_adn','pi'; 'mvdr_adn','mvdr'; 'mvdr_adn','pi_adn'; 'pi','single'; 'mvdr_adn','pi'};
for s = {'S1','S2','S3','all'}
  k = true(n,1); if ~strcmp(s{1},'all'), k = strcmp(sc,s{1}); end
  for q = 1:size(pairs,1)
    a = D.(pairs{q,1}).rmsJam(k); b = D.(pairs{q,2}).rmsJam(k); r = a./b; v = ~isnan(r);
    P('%-3s %-9s/%-7s n=%2d med %.3f [%.3f %.3f] better %2d/%2d\n', s{1}, pairs{q,1}, pairs{q,2}, sum(v), median(r(v)), prctile(r(v),25), prctile(r(v),75), sum(a(v)<b(v)), sum(v));
  end
end
P('\n== seed-3 ablation of the receiver aiding (ratio to pi)\n');
k = sd == 3;
for i = find(k)'
  P('%s_%s', sc{i}, pr{i}); for c = {'pi_a','pi_dn','pi_adn'}, P(' %s %.3f', c{1}, D.(c{1}).rmsJam(i)/D.pi.rmsJam(i)); end
  P(' | lost pi %.3f a %.3f dn %.3f adn %.3f\n', D.pi.lostFrac(i), D.pi_a.lostFrac(i), D.pi_dn.lostFrac(i), D.pi_adn.lostFrac(i));
end
for c = {'pi_a','pi_dn','pi_adn'}, P('median %s/pi %.3f, lost mean %.2f%%\n', c{1}, median(D.(c{1}).rmsJam(k)./D.pi.rmsJam(k),'omitnan'), 100*mean(D.(c{1}).lostFrac(k),'omitnan')); end
P('pi lost mean %.2f%%\n', 100*mean(D.pi.lostFrac(k)));
P('\n== max errors and coasting (UT)\n');
for c = main, P('%-9s maxJam med %.2f max %.2f | coastMax med %.2f (n %d)\n', c{1}, median(D.(c{1}).maxJam,'omitnan'), max(D.(c{1}).maxJam), median(D.(c{1}).coastMax,'omitnan'), sum(~isnan(D.(c{1}).coastMax))); end
clear cl
type(txt);

% ---- Fig. 9: UT rms during jamming per item, by scenario and configuration
rng(1);
col = {[0.32 0.32 0.31], [42 120 214]/255, [235 104 52]/255, [27 175 122]/255, [237 161 0]/255};
mk = {'x','o','s','^','d'};
fig = figure('Visible','off','Units','centimeters','Position',[0 0 17.5 6.2]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
scs = {'S1','S2','S3'};
for a = 1:3
  ax = nexttile(tl); hold(ax,'on'); set(ax,'YScale','log','FontSize',7,'Box','off','YGrid','on','GridAlpha',0.15,'TickDir','out');
  k = strcmp(sc, scs{a});
  for c = 1:5
    y = D.(main{c}).rmsJam(k);  x = c + 0.28*(rand(sum(k),1)-0.5);
    plot(ax, x, y, mk{c}, 'Color', col{c}, 'MarkerSize', 3.5, 'LineWidth', 0.6, 'MarkerFaceColor', 'none');
    m = median(y,'omitnan'); plot(ax, c + [-0.3 0.3], [m m], '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 1.2);
  end
  set(ax,'XLim',[0.4 5.6],'XTick',1:5,'XTickLabel',{'Single','PI','PI+aid','MVDR','MVDR+aid'},'YLim',[0.08 40]);
  xtickangle(ax, 30);
  title(ax, sprintf('%s (%d runs)', scs{a}, sum(k)), 'FontWeight','normal','FontSize',8);
  if a == 1, ylabel(ax, '3-D position error rms during jamming [m]'); end
end
set(fig,'PaperUnits','centimeters','PaperSize',[17.5 6.2],'PaperPosition',[0 0 17.5 6.2]);
print(fig, fullfile(outDir, 'fig_campaign.pdf'), '-dpdf', '-vector');
close(fig);
fprintf('paper_results: %s, %s\n', txt, fullfile(outDir, 'fig_campaign.pdf'));
end
