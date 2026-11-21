function run_campaign(worker, nWorkers, seeds)
%RUN_CAMPAIGN  Navigation-level scenarios S1-S3 (paper Table VIII), resumable.
%   run_campaign(worker, nWorkers, seeds)
%   Runs every (seed, scenario, profile) item whose index mod nWorkers equals
%   worker-1, in seed-major order, so that the first seeds of all scenarios
%   finish first. For each item: eight-element IF -> array processing in the
%   configurations of campaignConfig.m (single, pi, pi_adn, mvdr, mvdr_adn; S1 seeds 1-3
%   also pi_ff) -> receiver (NavIC,
%   CT + UT, warm start) -> metrics, saved to
%   J:\AJ-AS\code\campaign\results\<scen>_<profile>_s<seed>.mat.
%   Large IF files are deleted after each item. Items with results are skipped
%   (results lacking configurations added later are completed with finish_item);
%   items that fail are retried once after all the others.
%     matlab -batch "addpath('J:\AJ-AS\code\work\aj'); run_campaign(1,1,1:10)"
if nargin < 1, worker = 1; end
if nargin < 2, nWorkers = 1; end
if nargin < 3, seeds = 1:10; end
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'Simulator'), fullfile(repo, 'Simulator', 'F16Aircraft'), fullfile(repo, 'ArrayProcessor'), fullfile(repo, 'aj'));
root = ajPath('campaign');
prof = ajPath('data', 'FighterProfiles');
warmRef = ajPath('step4', 'clean_single', 'out');      % ephemeris decoded on clean data
resDir = fullfile(root, 'results');  if ~exist(resDir, 'dir'), mkdir(resDir); end
logf = fullfile(root, sprintf('campaign_worker%d.log', worker));

profiles = struct('S1', {{'baseline_cruise', 'rolls', 'break_gonset'}}, ...
                  'S2', {{'baseline_cruise', 'rolls'}}, 'S3', {{'baseline_cruise', 'rolls'}});
% profiles (trajectory + IMU), shared start with the baseline
for p = unique([profiles.S1 profiles.S2 profiles.S3])
    if ~exist(fullfile(prof, ['Fighter_' p{1} '_IMU_200Hz.mat']), 'file')
        here = pwd;  cd(fullfile(repo, 'Simulator'));
        generateFighterProfile(p{1}, struct('outDir', prof, 'accelRangeG', 20, 'seed', 1, 'T', 60));
        cd(here);
    end
end

items = {};
for s = seeds
    for sc = {'S1', 'S2', 'S3'}
        for p = profiles.(sc{1}), items(end+1, :) = {sc{1}, p{1}, s}; end %#ok<AGROW>
    end
end
mine = find(mod((1:size(items, 1)) - 1, nWorkers) == worker - 1);
logmsg(logf, sprintf('worker %d/%d: %d items', worker, nWorkers, numel(mine)));
tagRx = @(name) ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' name '_NavIC'];

for pass = 1:2              % second pass retries items that failed in the first
for i = mine
    [sc, name, seed] = items{i, :};
    id = sprintf('%s_%s_s%02d', sc, name, seed);
    % pick up code changed while this session runs (J: is a USB disk: MATLAB's
    % file-change notification does not see edits there)
    rehash path
    resFile = fullfile(resDir, [id '.mat']);
    [configs, spec] = campaignConfig(sc, seed);
    wroot = getenv('AJ_WORKDIR');  if isempty(wroot), wroot = fullfile(root, 'work'); end
    work = fullfile(wroot, id);
    if exist(resFile, 'file')
        % finished item: complete configurations added later (pi_adn, mvdr_adn, pi_dn)
        have = fieldnames(load(resFile));
        if all(ismember(configs, have)), continue; end
        logmsg(logf, ['complete ' id ' (' strjoin(setdiff(configs, have, 'stable'), ', ') ')']);
        try
            finish_item(id, work);
            logmsg(logf, ['completed ' id]);
        catch err
            logmsg(logf, sprintf('FAILED completing %s: %s (%s line %d)', id, err.message, err.stack(1).name, err.stack(1).line));
        end
        continue
    end
    tItem = tic;
    logmsg(logf, ['start ' id]);
    try
        [jam, info] = scenarioJammers(sc, seed);
        % temporary IF data (element files, processed streams): on a fast
        % disk when AJ_WORKDIR is set (e.g. the internal SSD); results stay on root
        g = generateIFdataArray(name, id, struct('profileDir', prof, 'outDir', fullfile(work, 'elements'), ...
            'ms', 60000, 'seed', 100 + seed, 'jammers', jam));
        popt = struct('profileDir', prof);
        if strcmp(sc, 'S3'), popt.blankK = 5; end
        R = struct('id', id, 'scen', sc, 'profile', name, 'seed', seed, 'info', info, 'gen', g.info);
        if strcmp(sc, 'S1'), jsFun = @(t) interp1([15 45], [40 100], min(max(t, 15), 45));
        else, jsFun = []; end
        for c = configs
            s = spec.(c{1});
            dd = fullfile(work, s.dataCfg, 'data');  oo = fullfile(work, c{1}, 'out');
            if ~exist(oo, 'dir'), mkdir(oo); end
            if isempty(s.aid)                % receiver-only configurations reuse dataCfg's streams
                ap = arrayProcess(g.files, g.truthFile, dd, name, s.amode, popt);
                R.(c{1}).arrayMinutes = ap.minutes;
            end
            for pre = {'eph', 'sbf', 'TckResultCT_forEph'}
                copyfile(fullfile(warmRef, [pre{1} '_' tagRx('baseline_cruise') '.mat']), ...
                         fullfile(oo, [pre{1} '_' tagRx(name) '.mat']));
            end
            RX = struct('dataDir', dd, 'outDir', oo, 'select', {{name}}, 'config', {{'NavIC'}}, 'mode', 3, ...
                        'rerun', true, 'arrayAid', s.aid);
            runReceiver(repo, RX);
            m = extractMetrics(oo, name, jsFun, struct('tOn', 15, 'tOff', 50));
            R.(c{1}).CT = m.CT;  R.(c{1}).UT = m.UT;
            R.(c{1}).arrayLog = decimateLog(load(fullfile(dd, ['array_log_' name '.mat'])));
            % keep the receiver results, delete the IF data after its last use
            if s.lastUse, rmdir(fullfile(dd, 'IFdata'), 's'); end
        end
        R.minutes = toc(tItem)/60;
        save(resFile, '-struct', 'R');
        rawIFStore('keep', id, fullfile(work, 'elements'));   % raw 8-element IF kept on J: for later re-processing
        % keep the receiver outputs (tables, plots, error CSVs) on the result disk
        if ~strncmpi(work, root, numel(root))
            keep = fullfile(root, 'work', id);  if ~exist(keep, 'dir'), mkdir(keep); end
            copyfile(work, keep);  rmdir(work, 's');
        end
        logmsg(logf, sprintf('done %s in %.1f min | UT rmsJam %s | UT tFail %s', id, R.minutes, ...
            strjoin(cellfun(@(c) sprintf('%s %.2f', c, R.(c).UT.rmsJam), configs, 'UniformOutput', false), ' '), ...
            strjoin(cellfun(@(c) sprintf('%s %g', c, R.(c).UT.tFail), configs, 'UniformOutput', false), ' ')));
    catch err
        logmsg(logf, sprintf('FAILED %s: %s (%s line %d)', id, err.message, err.stack(1).name, err.stack(1).line));
    end
end
end
logmsg(logf, 'worker finished');
end

function runReceiver(repo, RX) %#ok<INUSD>
% Receiver\init.m is a script that reads RX from the workspace
rehash path
run(fullfile(repo, 'Receiver', 'init.m'));
end

function logmsg(f, s)
fid = fopen(f, 'a');  fprintf(fid, '%s  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), s);  fclose(fid);
fprintf('%s\n', s);
end
