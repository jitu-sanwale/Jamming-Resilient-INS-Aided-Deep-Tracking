function export_localization(ids)
%EXPORT_LOCALIZATION Export truth and jammer data of stored raw-IF items for
%   aj/jammer_localization.py (the truth file holds a function handle that
%   Python cannot read).  Writes J:\AJ-AS\code\campaign\loc\loc_<id>.mat
%     export_localization({'S1_baseline_cruise_s01','S3_rolls_s02'})
%   Only stationary (ground) jammers are supported; S2 (airborne jammer) is skipped.
addpath(fileparts(mfilename('fullpath')));
outDir = ajPath('campaign', 'loc');  if ~exist(outDir, 'dir'), mkdir(outDir); end
for k = 1:numel(ids)
    id = ids{k};  tok = regexp(id, '^(S\d)_(.+)_s(\d+)$', 'tokens', 'once');
    if strcmp(tok{1}, 'S2'), fprintf('%s: airborne jammer, skipped\n', id); continue; end
    T = load(fullfile(ajPath('campaign', 'rawIF'), id, sprintf('navic8_truth_Fighter_%s_%s.mat', tok{2}, id)));
    jam = scenarioJammers(tok{1}, str2double(tok{3}));
    E = struct('t', T.t, 'att', T.att, 'posEcef', T.posEcef, 'uJamB', T.uJamB, 'uSatN', T.uSatN, 'prn', T.prn, ...
               'posB', T.array.posB, 'C', T.array.C, 'lambda', T.array.lambda);
    E.jamPos = reshape([jam.pos], 3, [])';  E.jamType = {jam.type};  E.jamOn = [jam.tOn];  E.jamOff = [jam.tOff];
    save(fullfile(outDir, ['loc_' id '.mat']), '-struct', 'E', '-v7');
    fprintf('%s: %d jammers %s\n', id, numel(jam), strjoin(E.jamType, ','));
end
end
