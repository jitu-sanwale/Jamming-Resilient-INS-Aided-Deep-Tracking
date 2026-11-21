function m = extractMetrics(outDir, name, jsFun, win)
%EXTRACTMETRICS Navigation metrics of one receiver run (CT and UT) for the result tables.
%   m = extractMetrics(outDir, name, jsFun, win)
%   outDir - receiver output folder of the run
%   name   - profile name
%   jsFun  - @(t) J/S [dB] of the strongest jammer at time t (S1 ramp), or []
%   win    - struct with tOn, tOff (jammer window, s), errThr (UT, m, default 10), errThrCT (CT, m, default 30),
%            truthFile (profile CSV; default ...\data\FighterProfiles\Fighter_<name>.csv)
%   m      - struct with fields CT and UT, each:
%              rmsPre, rmsJam, rmsPost, maxJam (3-D position error, m; windows
%              3.1 s-tOn, tOn-tOff, tOff+1 s-59.7 s), velRmsJam (m/s),
%              tFail (first time after tOn with error > threshold for 1 s or fewer
%              than four channels tracking; NaN if none), jsTol (J/S at tFail; Inf if
%              no failure), tRecover (time after tOff until the error stays
%              below the threshold and four channels track; NaN if not recovered).
%            CT: channels tracking = channels used by the CT EKF. CT errors are
%              computed here from navSolCT against the truth spline (the receiver
%              omits CT from its truth comparison in mode 3).
%            UT also: cn0Jam (mean C/N0 per channel during jamming), resets
%              (carrier-loop resets during jamming), chUsedMin (EKF-used channels),
%              chTrackMin / chTrackMean (channels with CaP reset in fewer than half the epochs of the last
%              1 s, during jamming), lostFrac (fraction of jamming epochs with
%              fewer than four tracking channels), coastMax (max 3-D error while
%              fewer than four channels track, i.e. INS bridging; NaN if never).
%            rxFailed: true when the receiver produced no solution (run aborted).
%   Failure uses the tracking-channel count: a carrier loop that is reset at
%   every epoch is not tracking even if its C/N0 estimate (noise floor ~25-35
%   dB-Hz) keeps it in the EKF.
if nargin < 4 || isempty(win), win = struct('tOn', 15, 'tOff', 50); end
if ~isfield(win, 'errThr'), win.errThr = 10; end          % UT
if ~isfield(win, 'errThrCT'), win.errThrCT = 30; end      % CT (clean CT error ~7 m rms)
if ~isfield(win, 'truthFile'), win.truthFile = fullfile(ajPath('data', 'FighterProfiles'), ['Fighter_' name '.csv']); end
rx0 = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Receiver');
if ~exist(rx0, 'dir'), rx0 = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Receiver'); end
addpath(fullfile(rx0, 'common'), fullfile(rx0, 'geo'));
res = fullfile(outDir, ['results_gps_0IF_Fs5MHz_8bit_IQ_Fighter_' name '_NavIC']);
tag = ['gps_0IF_Fs5MHz_8bit_IQ_Fighter_' name '_NavIC'];
m = struct();
for rx = {'CT', 'UT'}
    s = struct('rmsPre', NaN, 'rmsJam', NaN, 'rmsPost', NaN, 'maxJam', NaN, 'velRmsJam', NaN, ...
               'tFail', NaN, 'jsTol', NaN, 'tRecover', NaN);
    t = [];
    if strcmp(rx{1}, 'CT')
        f = fullfile(res, 'err_classical_10Hz.csv');
        nf = fullfile(outDir, ['navSolCT_' tag '.mat']);
        if exist(f, 'file')
            E = readmatrix(f);  t = E(:,1);  eP = E(:,2:4);  eV = E(:,5:7);
        elseif exist(nf, 'file') && exist(win.truthFile, 'file')
            N = load(nf, 'navSolutionsCT');  C = N.navSolutionsCT;
            tc = (100 + (1:size(C.usrPos, 1))')*1e-3;               % file.skip = 100 ms, as compareTruth
            i = 1:100:numel(tc);                                     % 10 Hz
            [eP, eV, ok] = errVsTruth(tc(i), C.usrPos(i,:), C.usrVel(i,:), win.truthFile);
            t = tc(i);  t = t(ok);  eP = eP(ok,:);  eV = eV(ok,:);
            nTrack = double(sum(C.usedCh(i(ok), :), 2));
        end
        if exist('nTrack', 'var'), trackOk = nTrack >= 4; else, trackOk = true(size(t)); end
    else
        f = fullfile(res, 'err_ultratight_10Hz.csv');
        nf = fullfile(outDir, ['navSolUT_' tag '.mat']);
        if exist(f, 'file')
            E = readmatrix(f);  t = E(:,1);  eP = E(:,2:4);  eV = E(:,5:7);
        elseif exist(nf, 'file') && exist(win.truthFile, 'file')
            % no truth comparison written (run stopped after tracking): errors
            % from the UT solution against the truth spline, at 10 Hz as compareTruth
            N = load(nf, 'navSolutionsUT');  U = N.navSolutionsUT;
            i = 1:max(1, round(0.1/median(diff(U.tData)))):numel(U.tData);
            tu = U.tData(i);  tu = tu(:);
            [eP, eV, ok] = errVsTruth(tu, U.usrPos(i,:), U.usrVel(i,:), win.truthFile);
            t = tu(ok);  eP = eP(ok,:);  eV = eV(ok,:);
        end
        trackOk = true(size(t));
        if ~isempty(t) && exist(nf, 'file')
            N = load(nf, 'navSolutionsUT');  U = N.navSolutionsUT;
            dt = median(diff(U.tData));  nw = max(1, round(1/dt));
            % tracking: CaP reset in fewer than half of the epochs of the last 1 s
            % (tracking loops: ~5 % of epochs; lost loops: 94-98 %)
            noReset = movmean(double(U.reset > 0), [nw-1 0], 1) < 0.5;
            nTr = sum(noReset & U.usedCh > 0, 2);
            trackOk = interp1(U.tData, double(nTr), t, 'previous', 'extrap') >= 4;
            j = U.tData >= win.tOn & U.tData < win.tOff;
            s.cn0Jam = mean(U.CN0(j, :), 1);
            s.resets = sum(U.reset(j, :), 'all');
            s.chUsedMin = min(sum(U.usedCh(j, :), 2));
            s.chTrackMin = min(nTr(j));
            s.chTrackMean = mean(nTr(j));
            s.lostFrac = mean(nTr(j) < 4);
        end
    end
    if ~isempty(t)
        e = sqrt(sum(eP.^2, 2));  v = sqrt(sum(eV.^2, 2));
        bad = ~isfinite(e);  e(bad) = 1e6;  v(~isfinite(v)) = 1e6;
        w = @(a, b) t >= a & t < b;
        s.rmsPre  = rmsOf(e(w(3.1, win.tOn)));
        s.rmsJam  = rmsOf(e(w(win.tOn, win.tOff)));
        s.rmsPost = rmsOf(e(w(win.tOff + 1, 59.7)));
        s.maxJam  = max([e(w(win.tOn, win.tOff)); NaN]);
        s.velRmsJam = rmsOf(v(w(win.tOn, win.tOff)));
        if strcmp(rx{1}, 'UT')
            c = w(win.tOn, win.tOff) & ~trackOk;
            s.coastMax = max([e(c); NaN]);
        end
        if strcmp(rx{1}, 'CT'), thr = win.errThrCT; else, thr = win.errThr; end
        over = e > thr;
        fail = (movmin(double(over), [0 9]) > 0) | ~trackOk;        % 1 s sustained (10 samples)
        k = find(fail & t >= win.tOn & t < win.tOff, 1);
        if ~isempty(k)
            s.tFail = t(k);
            if ~isempty(jsFun), s.jsTol = jsFun(t(k)); end
        elseif ~isempty(jsFun)
            s.jsTol = Inf;
        end
        kk = find(t >= win.tOff, 1);
        if ~isempty(kk)
            post = ~over(kk:end) & trackOk(kk:end);
            stay = flipud(cummin(flipud(post)));                    % stays good from here on
            r = find(stay, 1);
            if ~isempty(r), s.tRecover = t(kk + r - 1) - win.tOff; end
        end
    end
    s.rxFailed = isempty(t);        % receiver run produced no solution (aborted)
    m.(rx{1}) = s;
    clear nTrack
end
end

function [eP, eV, ok] = errVsTruth(t, pos, vel, truthFile)
% ENU position / velocity error against the profile truth spline (as compareTruth)
persistent cacheF cacheT
if isempty(cacheF) || ~strcmp(cacheF, truthFile)
    cacheT = trajectorySplineECEF(truthFile);  cacheF = truthFile;
end
traj = cacheT;
ok = t >= traj.t0 & t <= traj.t1;
eP = NaN(numel(t), 3);  eV = eP;
tt = t(ok);
rT = ppval(traj.ppPos, tt')';  vT = ppval(traj.ppVel, tt')';
dr = pos(ok,:) - rT;  dv = vel(ok,:) - vT;
llh = zeros(numel(tt), 2);
for i = 1:numel(tt), q = xyz2llh(rT(i,:)); llh(i,:) = q(1:2); end
sL = sin(llh(:,1)); cL = cos(llh(:,1)); sl = sin(llh(:,2)); cl = cos(llh(:,2));
rot = @(d) [-sl.*d(:,1) + cl.*d(:,2), ...
            -sL.*cl.*d(:,1) - sL.*sl.*d(:,2) + cL.*d(:,3), ...
             cL.*cl.*d(:,1) + cL.*sl.*d(:,2) + sL.*d(:,3)];
eP(ok,:) = rot(dr);  eV(ok,:) = rot(dv);
end

function r = rmsOf(x)
if isempty(x), r = NaN; else, r = sqrt(mean(x.^2)); end
end
