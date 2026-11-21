%RUN_STEP4_TESTS  First end-to-end runs: eight-element generator -> array processor -> receiver.
%   Clean seven-satellite record: single element (Tier 0 check against the
%   single-antenna seven-satellite result) and per-satellite MVDR beams.
%   Broadband ground jammer, J/S 60 dB from 10 s to 50 s: single element, PI
%   (shared beam) and per-satellite MVDR. Receiver: NavIC alone, CT + UT.
%     matlab -batch "run('J:\AJ-AS\code\work\aj\run_step4_tests.m')"
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'Simulator'), fullfile(repo, 'Simulator', 'F16Aircraft'), fullfile(repo, 'ArrayProcessor'));
prof = ajPath('data', 'FighterProfiles');
arr  = ajPath('arr');
root = ajPath('step4');
name = 'baseline_cruise';
tAll = tic;

%% generation
t = tic;
gc = generateIFdataArray(name, 'clean', struct('profileDir', prof, 'outDir', arr, 'ms', 60000, 'seed', 0));
jam = struct('type', 'bbn', 'pos', [19.85, 74.20, 0], 'jsDb', 60, 'tOn', 10, 'tOff', 50, 'params', struct());
gj = generateIFdataArray(name, 'bbn60', struct('profileDir', prof, 'outDir', arr, 'ms', 60000, 'seed', 0, 'jammers', jam));
fprintf('[timing] generation (2 records): %.1f min\n', toc(t)/60);

%% array processing and receiver runs
runs = {gc, 'clean_single', 'single';  gc, 'clean_mvdr', 'mvdr'; ...
        gj, 'bbn60_single', 'single';  gj, 'bbn60_pi', 'pi';  gj, 'bbn60_mvdr', 'mvdr'};
tagRx = ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' name '_NavIC'];
refOut = fullfile(root, 'clean_single', 'out');
for k = 1:size(runs, 1)
    g = runs{k, 1};  id = runs{k, 2};  mode = runs{k, 3};
    dd = fullfile(root, id, 'data');  oo = fullfile(root, id, 'out');
    if ~exist(oo, 'dir'), mkdir(oo); end
    t = tic;
    arrayProcess(g.files, g.truthFile, dd, name, mode, struct('profileDir', prof));
    fprintf('[timing] array processing %s: %.1f min\n', id, toc(t)/60);
    if k > 1   % warm start: ephemeris decoded in the clean single-element run
        for pre = {'eph', 'sbf', 'TckResultCT_forEph'}
            copyfile(fullfile(refOut, [pre{1} '_' tagRx '.mat']), fullfile(oo, [pre{1} '_' tagRx '.mat']));
        end
    end
    RX = struct('dataDir', dd, 'outDir', oo, 'select', {{name}}, 'config', {{'NavIC'}}, 'mode', 3, 'rerun', true);
    t = tic;
    run(fullfile(repo, 'Receiver', 'init.m'));
    fprintf('[timing] receiver %s: %.1f min\n', id, toc(t)/60);
end
fprintf('[timing] total: %.1f min\n', toc(tAll)/60);
