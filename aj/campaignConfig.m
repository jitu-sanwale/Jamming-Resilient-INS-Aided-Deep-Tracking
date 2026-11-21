function [configs, spec] = campaignConfig(sc, seed)
%CAMPAIGNCONFIG Receiver/array configurations of one campaign item.
%   [configs, spec] = campaignConfig(sc, seed)
%   configs - cell array in processing order
%   spec    - struct per configuration: amode (arrayProcess mode), aid (receiver
%             CRPA aiding, initParameters track.utArrayAid), dataCfg (configuration
%             whose processed streams it reads), lastUse (true for the last
%             configuration that reads dataCfg's streams, which are then deleted)
%     single    element 0 (single-antenna reference)
%     pi        power inversion, shared beam (paper C1 with CT, C2 with UT)
%     pi_adn    pi streams + receiver carrier feed-forward of the beam phase (a),
%               array-informed C/N0 cap with an envelope-tracked reference (d) and
%               blanking-aware loop C/N0 (n)
%     pi_a, pi_dn  ablation on seed 3: (a) only, (d)+(n) only
%     mvdr      INS-steered per-satellite MVDR beams (paper C3)
%     mvdr_adn  mvdr streams + (a), (d) and (n)
%     pi_ff     PI with per-satellite phase compensation in the array (S1, seeds 1-3)
%   pi_ab, mvdr_ab, pi_b (the first weighting, 'b': cap calibrated over the first
%   2 s of UT) remain in the results of the items run before 2026-11-08; that
%   cap is biased low when interference is present during calibration (S3).
configs = {'single', 'pi', 'pi_adn', 'mvdr', 'mvdr_adn'};
if seed == 3, configs = [configs(1:3), {'pi_a', 'pi_dn'}, configs(4:end)]; end   % ablation of (a) and (d)+(n)
if strcmp(sc, 'S1') && seed <= 3, configs{end+1} = 'pi_ff'; end
spec = struct();
for c = configs
    k = c{1};
    tok = regexp(k, '^(pi|mvdr)_(adn|dn|ab|a|b)$', 'tokens', 'once');      % receiver-aided variants
    if ~isempty(tok), s = struct('amode', tok{1}, 'aid', tok{2});
    else, s = struct('amode', k, 'aid', ''); end
    s.dataCfg = s.amode;
    spec.(k) = s;
end
for c = configs
    users = configs(cellfun(@(x) strcmp(spec.(x).dataCfg, spec.(c{1}).dataCfg), configs));
    spec.(c{1}).lastUse = strcmp(users{end}, c{1});
end
end
