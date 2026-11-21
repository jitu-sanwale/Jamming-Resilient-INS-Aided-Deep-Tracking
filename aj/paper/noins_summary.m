function T = noins_summary(outCsv)
%NOINS_SUMMARY  INS contribution to the beamformer (R1 Section VII-D).
%   T = noins_summary(outCsv)
%   Compares, run by run, MVDR steered without INS (mvdr_ni: attitude from the
%   GNSS velocity, roll 0, fixed 10-ms block; run_noins_campaign) with the
%   INS-steered MVDR (mvdr) and PI of run_campaign, from
%   ajPath('campaign','results'): noins_<id>.mat and <id>.mat. Prints the group
%   statistics quoted in Section VII-D and writes the per-run table to outCsv
%   (default: ajPath('campaign','noins_summary.csv')).
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo, 'aj'));
R = ajPath('campaign', 'results');
if nargin < 1, outCsv = ajPath('campaign', 'noins_summary.csv'); end
d = dir(fullfile(R, 'noins_S*.mat'));
rows = {};
for i = 1:numel(d)
    A = load(fullfile(R, d(i).name));  id = A.id;  B = load(fullfile(R, [id '.mat']));
    tok = regexp(id, '^(S\d)_(.+)_s(\d+)$', 'tokens', 'once');
    a = A.mvdr_ni.UT;  b = B.mvdr.UT;  p = B.pi.UT;
    L = A.mvdr_ni.arrayLog;  w = L.t > 15 & L.t < 50;
    gl = NaN;  if isfield(L, 'gainTrueDb'), gl = median(min(L.gainTrueDb(w, :), [], 1)); end   % median over satellites of the worst true-direction gain
    rows(end+1, :) = {tok{1}, tok{2}, str2double(tok{3}), a.rmsJam, b.rmsJam, p.rmsJam, a.lostFrac, b.lostFrac, ...
        ~isnan(a.tFail), ~isnan(b.tFail), mean(a.cn0Jam), mean(b.cn0Jam), a.resets, b.resets, a.chTrackMean, b.chTrackMean, gl}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'scen','prof','seed','rmsNI','rmsINS','rmsPI','lostNI','lostINS','failNI','failINS','cn0NI','cn0INS','resNI','resINS','chNI','chINS','worstGainTrueDb'});
T.man = ~strcmp(T.prof, 'baseline_cruise');
fprintf('runs: %d (cruise %d, maneuver %d)\n', height(T), sum(~T.man), sum(T.man));
grp = {'cruise', ~T.man; 'maneuver', T.man; 'S1', strcmp(T.scen,'S1'); 'S2', strcmp(T.scen,'S2'); 'S3', strcmp(T.scen,'S3'); ...
       'S1 rolls', strcmp(T.scen,'S1')&strcmp(T.prof,'rolls'); 'S1 break', strcmp(T.prof,'break_gonset'); ...
       'S2 rolls', strcmp(T.scen,'S2')&strcmp(T.prof,'rolls'); 'S3 rolls', strcmp(T.scen,'S3')&strcmp(T.prof,'rolls'); 'all', true(height(T),1)};
fprintf('%-9s %3s | rms med NI/INS/PI | p95 NI/INS | ratio NI/INS (NI worse) | lost%% NI/INS | runs lost>0 NI/INS | fail NI/INS | ch NI/INS | C/N0 NI-INS | resets NI/INS | worst true gain dB\n', 'group', 'n');
for g = 1:size(grp, 1)
    k = grp{g, 2};  t = T(k, :);  r = t.rmsNI ./ t.rmsINS;
    fprintf('%-9s %3d | %.2f/%.2f/%.2f | %.2f/%.2f | %.2f (%d/%d) | %.2f/%.2f | %d/%d | %d/%d | %.2f/%.2f | %+.1f | %.0f/%.0f | %.1f\n', grp{g,1}, height(t), ...
        median(t.rmsNI), median(t.rmsINS), median(t.rmsPI), prctile(t.rmsNI,95), prctile(t.rmsINS,95), median(r), sum(r > 1), height(t), ...
        100*mean(t.lostNI), 100*mean(t.lostINS), sum(t.lostNI > 0), sum(t.lostINS > 0), sum(t.failNI), sum(t.failINS), ...
        mean(t.chNI), mean(t.chINS), mean(t.cn0NI - t.cn0INS), mean(t.resNI), mean(t.resINS), median(t.worstGainTrueDb, 'omitnan'));
end
fprintf('maneuver lostNI range %.1f-%.1f%%, lostINS max %.1f%%; NI worse than PI in rms: %d runs\n', 100*min(T.lostNI(T.man)), 100*max(T.lostNI(T.man)), 100*max(T.lostINS), sum(T.rmsNI > T.rmsPI));
% two-sided sign test on the paired time with fewer than four channels (maneuvers)
k = T.man; n = sum(T.lostNI(k) ~= T.lostINS(k)); x = sum(T.lostNI(k) > T.lostINS(k));
pval = 2*sum(arrayfun(@(j) nchoosek(n, j), x:n))/2^n;  fprintf('sign test maneuvers lostNI>lostINS: %d of %d, p = %.2g\n', x, n, min(pval,1));
writetable(T, outCsv);
end
