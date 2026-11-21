function T = insfix_summary(outCsv)
%INSFIX_SUMMARY  Separates INS steering from the body-rate schedule (R1 Section VII-D).
%   T = insfix_summary(outCsv)
%   Compares, run by run, three MVDR configurations from ajPath('campaign','results'):
%     IF  = INS steering + fixed 10-ms block (mvdr_if,  insfix_<id>.mat; run_insfix_campaign)
%     INS = INS steering + body-rate schedule (mvdr,    <id>.mat;        run_campaign)
%     NI  = GNSS steering + fixed 10-ms block (mvdr_ni, noins_<id>.mat;  run_noins_campaign)
%   IF vs INS isolates the schedule; NI vs IF isolates the steering. Prints the
%   group statistics and writes the per-run table to outCsv
%   (default: ajPath('campaign','insfix_summary.csv')).
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo, 'aj'));
R = ajPath('campaign', 'results');
if nargin < 1, outCsv = ajPath('campaign', 'insfix_summary.csv'); end
d = dir(fullfile(R, 'insfix_S*.mat'));
rows = {};
for i = 1:numel(d)
    F = load(fullfile(R, d(i).name));  id = F.id;
    B = load(fullfile(R, [id '.mat']));  N = load(fullfile(R, ['noins_' id '.mat']));
    tok = regexp(id, '^(S\d)_(.+)_s(\d+)$', 'tokens', 'once');
    f = F.mvdr_if.UT;  b = B.mvdr.UT;  a = N.mvdr_ni.UT;
    rows(end+1, :) = {tok{1}, tok{2}, str2double(tok{3}), ...
        f.rmsJam, b.rmsJam, a.rmsJam, f.lostFrac, b.lostFrac, a.lostFrac, ...
        mean(f.cn0Jam), mean(b.cn0Jam), mean(a.cn0Jam), f.resets, b.resets, a.resets, ...
        f.chTrackMean, b.chTrackMean, a.chTrackMean}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'scen','prof','seed', 'rmsIF','rmsINS','rmsNI', 'lostIF','lostINS','lostNI', ...
    'cn0IF','cn0INS','cn0NI', 'resIF','resINS','resNI', 'chIF','chINS','chNI'});
T.man = ~strcmp(T.prof, 'baseline_cruise');
fprintf('runs: %d (cruise %d, maneuver %d)\n', height(T), sum(~T.man), sum(T.man));
grp = {'cruise', ~T.man; 'maneuver', T.man; 'S1 rolls', strcmp(T.scen,'S1')&strcmp(T.prof,'rolls'); ...
       'S1 break', strcmp(T.prof,'break_gonset'); 'S2 rolls', strcmp(T.scen,'S2')&strcmp(T.prof,'rolls'); ...
       'S3 rolls', strcmp(T.scen,'S3')&strcmp(T.prof,'rolls')};
fprintf('%-9s %3s | rms med IF/INS/NI | lost%% IF/INS/NI | runs lost>0 IF/INS/NI | ch IF/INS/NI | resets IF/INS/NI | C/N0 IF-INS, NI-IF | rms ratio IF/INS (IF worse), NI/IF (NI worse)\n', 'group', 'n');
for g = 1:size(grp, 1)
    t = T(grp{g, 2}, :);  r1 = t.rmsIF ./ t.rmsINS;  r2 = t.rmsNI ./ t.rmsIF;
    fprintf('%-9s %3d | %.2f/%.2f/%.2f | %.2f/%.2f/%.2f | %d/%d/%d | %.2f/%.2f/%.2f | %.1f/%.1f/%.1f | %+.2f, %+.2f | %.2f (%d/%d), %.2f (%d/%d)\n', ...
        grp{g,1}, height(t), median(t.rmsIF), median(t.rmsINS), median(t.rmsNI), ...
        100*mean(t.lostIF), 100*mean(t.lostINS), 100*mean(t.lostNI), sum(t.lostIF>0), sum(t.lostINS>0), sum(t.lostNI>0), ...
        mean(t.chIF), mean(t.chINS), mean(t.chNI), mean(t.resIF), mean(t.resINS), mean(t.resNI), ...
        mean(t.cn0IF - t.cn0INS), mean(t.cn0NI - t.cn0IF), median(r1), sum(r1 > 1), height(t), median(r2), sum(r2 > 1), height(t));
end
k = T.man;
fprintf('maneuver lost%% range IF %.2f-%.2f, INS %.2f-%.2f, NI %.1f-%.1f\n', 100*min(T.lostIF(k)), 100*max(T.lostIF(k)), ...
    100*min(T.lostINS(k)), 100*max(T.lostINS(k)), 100*min(T.lostNI(k)), 100*max(T.lostNI(k)));
signtest2('lostIF>lostINS (schedule)', T.lostIF(k), T.lostINS(k));
signtest2('lostNI>lostIF (steering)',  T.lostNI(k), T.lostIF(k));
signtest2('rmsIF>rmsINS (schedule)',   T.rmsIF(k),  T.rmsINS(k));
signtest2('rmsNI>rmsIF (steering)',    T.rmsNI(k),  T.rmsIF(k));
writetable(T, outCsv);
end

function signtest2(name, x, y)
% two-sided exact sign test on paired values, ties dropped
n = sum(x ~= y);  m = sum(x > y);
if n == 0, fprintf('sign test maneuvers %s: no non-tied pairs\n', name); return; end
p = 2*sum(arrayfun(@(j) nchoosek(n, j), max(m, n-m):n))/2^n;
fprintf('sign test maneuvers %s: %d of %d non-tied, p = %.2g\n', name, m, n, min(p, 1));
end
