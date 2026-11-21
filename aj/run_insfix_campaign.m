function run_insfix_campaign(worker, nWorkers, seeds)
%RUN_INSFIX_CAMPAIGN  INS-steered MVDR with a fixed 10-ms block (mvdr_if), resumable.
%   run_insfix_campaign(worker, nWorkers, seeds)
%   Separates the two effects of the no-INS ablation (run_noins_campaign): the
%   beams are steered with the INS attitude, exactly as in run_campaign's 'mvdr'
%   configuration, but the covariance block is fixed at 10 ms instead of the
%   body-rate schedule (8). Compared with 'mvdr' (INS steering + schedule) it
%   isolates the scheduling; compared with 'mvdr_ni' (GNSS-velocity steering +
%   10 ms) it isolates the steering. Same items and worker split as run_campaign.
%   Each item: stored 8-element IF (ajPath('campaign','rawIF',<id>)) -> arrayProcess
%   'mvdr' with tbMs = 10 -> deep-tracking receiver without array aiding ->
%   results\insfix_<id>.mat. Items with a result are skipped.
%     matlab -batch "addpath('J:\AJ-AS\code\work\aj'); run_insfix_campaign(1,2,1:10)"
if nargin < 1, worker = 1; end
if nargin < 2, nWorkers = 1; end
if nargin < 3, seeds = 1:10; end
root = ajPath('campaign');
logf = fullfile(root, sprintf('insfix_worker%d.log', worker));
profiles = struct('S1', {{'baseline_cruise', 'rolls', 'break_gonset'}}, ...
                  'S2', {{'baseline_cruise', 'rolls'}}, 'S3', {{'baseline_cruise', 'rolls'}});
items = {};
for s = seeds
    for sc = {'S1', 'S2', 'S3'}
        for p = profiles.(sc{1}), items{end+1} = sprintf('%s_%s_s%02d', sc{1}, p{1}, s); end %#ok<AGROW>
    end
end
mine = find(mod((1:numel(items)) - 1, nWorkers) == worker - 1);
logmsg(logf, sprintf('worker %d/%d: %d items', worker, nWorkers, numel(mine)));
for pass = 1:2              % second pass retries failed items
    for i = mine
        id = items{i};
        if exist(fullfile(root, 'results', ['insfix_' id '.mat']), 'file'), continue; end
        rehash path
        logmsg(logf, ['start ' id]);
        try
            line = insfixItem(id, 10);
            logmsg(logf, line);
        catch err
            logmsg(logf, sprintf('FAILED %s: %s (%s line %d)', id, err.message, err.stack(1).name, err.stack(1).line));
        end
    end
end
logmsg(logf, 'worker finished');
end

function line = insfixItem(id, tbMs)
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'Simulator'), fullfile(repo, 'ArrayProcessor'), fullfile(repo, 'aj'));
root = ajPath('campaign');  prof = ajPath('data', 'FighterProfiles');
warmRef = ajPath('step4', 'clean_single', 'out');
tok = regexp(id, '^(S\d)_(.+)_s(\d+)$', 'tokens', 'once');
sc = tok{1};  name = tok{2};
el = fullfile(root, 'rawIF', id);
files = arrayfun(@(m) fullfile(el, sprintf('navic8_e%d_0IF_Fs5MHz_16bit_IQ_Fighter_%s_%s.bin', m, name, id)), 0:7, 'UniformOutput', false);
truthFile = fullfile(el, sprintf('navic8_truth_Fighter_%s_%s.mat', name, id));
assert(exist(truthFile, 'file') == 2, 'no stored raw IF for %s', id);
wroot = getenv('AJ_WORKDIR');  if isempty(wroot), wroot = fullfile(root, 'insfix'); end
work = fullfile(wroot, id);  dd = fullfile(work, 'mvdr_if', 'data');  oo = fullfile(work, 'mvdr_if', 'out');
if ~exist(oo, 'dir'), mkdir(oo); end
tagRx = @(n) ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' n '_NavIC'];
tItem = tic;
popt = struct('profileDir', prof, 'steerAtt', 'ins', 'tbMs', tbMs);   % INS steering (as 'mvdr'), fixed block
if strcmp(sc, 'S3'), popt.blankK = 5; end
ap = arrayProcess(files, truthFile, dd, name, 'mvdr', popt);
for pre = {'eph', 'sbf', 'TckResultCT_forEph'}
    copyfile(fullfile(warmRef, [pre{1} '_' tagRx('baseline_cruise') '.mat']), fullfile(oo, [pre{1} '_' tagRx(name) '.mat']));
end
runReceiver(repo, struct('dataDir', dd, 'outDir', oo, 'select', {{name}}, 'config', {{'NavIC'}}, 'mode', 3, 'rerun', true, 'arrayAid', ''));
if strcmp(sc, 'S1'), jsFun = @(t) interp1([15 45], [40 100], min(max(t, 15), 45)); else, jsFun = []; end
m = extractMetrics(oo, name, jsFun, struct('tOn', 15, 'tOff', 50));
R = struct('id', id, 'tbMs', tbMs, 'steerAtt', 'ins', 'arrayMinutes', ap.minutes);
R.mvdr_if.CT = m.CT;  R.mvdr_if.UT = m.UT;
R.mvdr_if.arrayLog = decimateLog(load(fullfile(dd, ['array_log_' name '.mat'])));
R.minutes = toc(tItem)/60;
save(fullfile(root, 'results', ['insfix_' id '.mat']), '-struct', 'R');
% keep receiver outputs on J:, drop the processed streams
if exist(fullfile(dd, 'IFdata'), 'dir'), rmdir(fullfile(dd, 'IFdata'), 's'); end
keep = fullfile(root, 'insfix', id);
if ~strcmpi(work, keep)
    if ~exist(keep, 'dir'), mkdir(keep); end
    copyfile(work, keep);  rmdir(work, 's');
end
B = load(fullfile(root, 'results', [id '.mat']));
N = load(fullfile(root, 'results', ['noins_' id '.mat']));
line = sprintf('done %s in %.1f min | UT rmsJam if %.2f ins %.2f ni %.2f | lostFrac if %.3f ins %.3f ni %.3f', ...
    id, R.minutes, m.UT.rmsJam, B.mvdr.UT.rmsJam, N.mvdr_ni.UT.rmsJam, m.UT.lostFrac, B.mvdr.UT.lostFrac, N.mvdr_ni.UT.lostFrac);
end

function runReceiver(repo, RX) %#ok<INUSD>
rehash path
run(fullfile(repo, 'Receiver', 'init.m'));
end

function logmsg(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), s);  fclose(fid);
fprintf('%s\n', s);
end
