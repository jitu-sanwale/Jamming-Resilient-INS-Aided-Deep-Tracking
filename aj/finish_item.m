function finish_item(id, work)
%FINISH_ITEM  Complete a campaign item that failed part-way, reusing what is on disk.
%   finish_item('S1_rolls_s01', 'C:\ajtmp\S1_rolls_s01')
%   Configurations whose receiver output exists (navSolUT_*.mat in work\<cfg>\out
%   or in the kept folder J:\AJ-AS\code\campaign\work\<id>\<cfg>\out) are only
%   re-scored; the others are array-processed from work\elements (fetched from
%   the raw-IF store J:\AJ-AS\code\campaign\rawIF, or regenerated with the
%   item's seed when absent, i.e. identical data) and run through the
%   receiver. Then the result file is written (fields of an existing one kept)
%   and the work folder merged into the kept folder, as run_campaign does.
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'Simulator'), fullfile(repo, 'Simulator', 'F16Aircraft'), fullfile(repo, 'ArrayProcessor'), fullfile(repo, 'aj'));
root = ajPath('campaign');  prof = ajPath('data', 'FighterProfiles');
warmRef = ajPath('step4', 'clean_single', 'out');
tok = regexp(id, '^(S\d)_(.+)_s(\d+)$', 'tokens', 'once');
sc = tok{1};  name = tok{2};  seed = str2double(tok{3});
[~, info] = scenarioJammers(sc, seed);
tagRx = @(n) ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' n '_NavIC'];
el = fullfile(work, 'elements');
keep = fullfile(root, 'work', id);
needGen = ~exist(fullfile(el, sprintf('navic8_truth_Fighter_%s_%s.mat', name, id)), 'file');
files = arrayfun(@(m) fullfile(el, sprintf('navic8_e%d_0IF_Fs5MHz_16bit_IQ_Fighter_%s_%s.bin', m, name, id)), 0:7, 'UniformOutput', false);
truthFile = fullfile(el, sprintf('navic8_truth_Fighter_%s_%s.mat', name, id));
[configs, spec] = campaignConfig(sc, seed);
popt = struct('profileDir', prof);  if strcmp(sc, 'S3'), popt.blankK = 5; end
if strcmp(sc, 'S1'), jsFun = @(t) interp1([15 45], [40 100], min(max(t, 15), 45)); else, jsFun = []; end
resFile = fullfile(root, 'results', [id '.mat']);
if exist(resFile, 'file'), R = load(resFile);
else, R = struct('id', id, 'scen', sc, 'profile', name, 'seed', seed, 'info', info, 'gen', struct('note', 'finished by finish_item'));
end
tItem = tic;
streamOf = @(d) fullfile(d, 'IFdata', ['navic_0IF_Fs5MHz_8bit_IQ_Fighter_' name '.bin']);
for c = configs
    s = spec.(c{1});
    dd = fullfile(work, s.dataCfg, 'data');  oo = fullfile(work, c{1}, 'out');
    done = exist(fullfile(oo, ['navSolUT_' tagRx(name) '.mat']), 'file') == 2;
    if ~done && exist(fullfile(keep, c{1}, 'out', ['navSolUT_' tagRx(name) '.mat']), 'file')
        dd = fullfile(keep, s.dataCfg, 'data');  oo = fullfile(keep, c{1}, 'out');  done = true;
    end
    if ~done
        if ~exist(streamOf(dd), 'file')
            % streams of dataCfg: process them (from regenerated elements if needed)
            if needGen && rawIFStore('fetch', id, el)          % stored raw IF instead of regenerating it
                needGen = false;
            end
            if needGen
                jam = scenarioJammers(sc, seed);
                g = generateIFdataArray(name, id, struct('profileDir', prof, 'outDir', el, 'ms', 60000, 'seed', 100 + seed, 'jammers', jam));
                if ~isfield(R, 'gen') || isfield(R.gen, 'note'), R.gen = g.info; end
                needGen = false;
            end
            ap = arrayProcess(files, truthFile, dd, name, s.amode, popt);
            if isempty(s.aid), R.(c{1}).arrayMinutes = ap.minutes; end
        end
        if ~exist(oo, 'dir'), mkdir(oo); end
        for pre = {'eph', 'sbf', 'TckResultCT_forEph'}
            copyfile(fullfile(warmRef, [pre{1} '_' tagRx('baseline_cruise') '.mat']), fullfile(oo, [pre{1} '_' tagRx(name) '.mat']));
        end
        runReceiver(repo, struct('dataDir', dd, 'outDir', oo, 'select', {{name}}, 'config', {{'NavIC'}}, 'mode', 3, ...
                                 'rerun', true, 'arrayAid', s.aid));
    end
    m = extractMetrics(oo, name, jsFun, struct('tOn', 15, 'tOff', 50));
    R.(c{1}).CT = m.CT;  R.(c{1}).UT = m.UT;
    R.(c{1}).arrayLog = decimateLog(load(fullfile(dd, ['array_log_' name '.mat'])));
    if s.lastUse && exist(fullfile(dd, 'IFdata'), 'dir'), rmdir(fullfile(dd, 'IFdata'), 's'); end
end
R.minutes = toc(tItem)/60;  R.metricsVersion = 2;
save(resFile, '-struct', 'R');
if exist(el, 'dir')
    if rawIFStore('has', id), rmdir(el, 's'); else, rawIFStore('keep', id, el); end
end
if ~exist(keep, 'dir'), mkdir(keep); end
if ~strcmpi(work, keep), copyfile(work, keep);  rmdir(work, 's'); end
fprintf('done %s | UT rmsJam %s\n', id, strjoin(cellfun(@(c) sprintf('%s %.2f', c, R.(c).UT.rmsJam), configs, 'UniformOutput', false), ', '));
end

function runReceiver(repo, RX) %#ok<INUSD>
% Receiver\init.m is a script that reads RX from the workspace
rehash path
run(fullfile(repo, 'Receiver', 'init.m'));
end
