function ajMenu(AJ)
%AJMENU  Mode selection for the CRPA anti-jamming extension (IEEE TAES R1).
%   ajMenu        asks for the mode and its options in the command window.
%   ajMenu(AJ)    runs without prompts; fields of AJ (all optional):
%     mode      1-8, see the menu below (default: ask)
%     worker    worker index for the campaign modes (default 1)
%     nWorkers  number of parallel MATLAB sessions (default 1)
%     seeds     realizations, e.g. 1:10 (default 1:10)
%     ids       item ids for jammer localization, e.g.
%               {'S1_baseline_cruise_s01', 'S1_break_gonset_s01'}
%
%   Modes (run them in this order to reproduce the manuscript; see README):
%     1) jammer-free reference and first jammed checks   (run_step4_tests)
%     2) navigation-level campaign, S1-S3                (run_campaign)
%     3) no-INS ablation: GNSS steering, 10-ms block     (run_noins_campaign)
%     4) INS steering with a fixed 10-ms block           (run_insfix_campaign)
%     5) tables, Fig. 9 and Section VII-D statistics     (paper_results, noins_summary, insfix_summary)
%     6) jammer localization, Fig. 8                     (export_localization, jammer_localization.py)
%     7) array-level figures, Figs. 3, 6, 7 and 8        (Python scripts in aj/paper/figures)
%     8) one realization: modes 2-4 for a single seed (the 7 items S1-S3 of that seed)
%
%   The campaign modes are resumable: items that already have a result are
%   skipped. Data folders: see ajPath (environment variable AJ_DATA_ROOT).
%
%   Example, worker 2 of 2 of the main campaign without prompts:
%     ajMenu(struct('mode', 2, 'worker', 2, 'nWorkers', 2, 'seeds', 1:10))

if nargin < 1 || ~isstruct(AJ), AJ = struct(); end
here = fileparts(mfilename('fullpath'));          % <repo>/aj
repo = fileparts(here);
addpath(here, fullfile(here, 'paper'));
ask = ~isfield(AJ, 'mode');

modes = {'jammer-free reference and first jammed checks (run_step4_tests)'
         'navigation-level campaign, scenarios S1-S3 (run_campaign)'
         'no-INS ablation: GNSS steering, fixed 10-ms block (run_noins_campaign)'
         'INS steering with a fixed 10-ms block (run_insfix_campaign)'
         'tables, Fig. 9 and Section VII-D statistics (from campaign results)'
         'jammer localization, Fig. 8 (needs stored raw IF and Python 3)'
         'array-level figures, Figs. 3, 6, 7 and 8 (Python 3)'
         'one realization: modes 2-4 for a single seed (7 items)'};
if ask
    fprintf('\nCRPA anti-jamming, NavIC L5, eight-element array (IEEE TAES R1)\n');
    for k = 1:numel(modes), fprintf('  %d) %s\n', k, modes{k}); end
    AJ.mode = str2double(strtrim(input('Select mode: ', 's')));
end
if ~isscalar(AJ.mode) || ~any(AJ.mode == 1:numel(modes))
    error('ajMenu: choose a mode between 1 and %d.', numel(modes));
end

% options of the campaign modes
if any(AJ.mode == [2 3 4 8])
    if AJ.mode == 8
        AJ.worker = 1;  AJ.nWorkers = 1;
        AJ.seeds = getOpt(AJ, 'seeds', 1, ask, 'Seed (realization 1-10) [1]: ');
        AJ.seeds = AJ.seeds(1);
    else
        AJ.worker   = getOpt(AJ, 'worker', 1, ask, 'Worker index [1]: ');
        AJ.nWorkers = getOpt(AJ, 'nWorkers', 1, ask, 'Number of workers [1]: ');
        AJ.seeds    = getOpt(AJ, 'seeds', 1:10, ask, 'Seeds, e.g. 1:10 [1:10]: ');
    end
    if AJ.worker < 1 || AJ.worker > AJ.nWorkers
        error('ajMenu: worker must lie between 1 and nWorkers.');
    end
end

switch AJ.mode
    case 1
        run(fullfile(here, 'run_step4_tests.m'));
    case 2
        run_campaign(AJ.worker, AJ.nWorkers, AJ.seeds);
    case 3
        run_noins_campaign(AJ.worker, AJ.nWorkers, AJ.seeds);
    case 4
        run_insfix_campaign(AJ.worker, AJ.nWorkers, AJ.seeds);
    case 5
        paper_results;  noins_summary;  insfix_summary;
    case 6
        def = {'S1_baseline_cruise_s01', 'S1_break_gonset_s01'};
        if isfield(AJ, 'ids'), ids = AJ.ids;
        elseif ask
            s = strtrim(input(sprintf('Item ids, separated by spaces [%s]: ', strjoin(def, ' ')), 's'));
            if isempty(s), ids = def; else, ids = strsplit(s); end
        else, ids = def;
        end
        export_localization(ids);
        runPython(repo, [{'aj/jammer_localization.py'} ids(:)']);
        jammerLocFig(repo);
    case 7
        for f = {'crpa8_element_patterns.py', 'crpa8_pattern.py', 'crpa8_roll_nulling.py'}
            runPython(repo, {['aj/paper/figures/' f{1}]});
        end
        jammerLocFig(repo);
    case 8
        % the three campaign drivers share the item list; one seed, one worker
        run_campaign(1, 1, AJ.seeds);
        run_noins_campaign(1, 1, AJ.seeds);
        run_insfix_campaign(1, 1, AJ.seeds);
end
end

function v = getOpt(AJ, name, def, ask, prompt)
% field of AJ, else the answer to a prompt, else the default
if isfield(AJ, name), v = AJ.(name); return; end
v = def;
if ask
    s = strtrim(input(prompt, 's'));
    if ~isempty(s), v = eval(s); end
end
end

function jammerLocFig(repo)
% Fig. 8 (fig_jammer_loc.py plots the two S1 items of the default list):
% two-column version, then the one-column version used in R1 (fig_jammer_loc_col)
old = getenv('FIG_LAYOUT');  c = onCleanup(@() setenv('FIG_LAYOUT', old));
setenv('FIG_LAYOUT', '');        runPython(repo, {'aj/paper/figures/fig_jammer_loc.py'});
setenv('FIG_LAYOUT', 'column');  runPython(repo, {'aj/paper/figures/fig_jammer_loc.py'});
end

function runPython(repo, args)
% run a Python 3 script from the repository folder; stop on failure
py = getenv('AJ_PYTHON');  if isempty(py), py = 'python'; end
cmd = sprintf('"%s" %s', py, strjoin(cellfun(@(a) ['"' a '"'], args, 'UniformOutput', false), ' '));
here = pwd;  cd(repo);  c = onCleanup(@() cd(here));
fprintf('%s\n', cmd);
if system(cmd) ~= 0, error('ajMenu: %s failed.', args{1}); end
end
