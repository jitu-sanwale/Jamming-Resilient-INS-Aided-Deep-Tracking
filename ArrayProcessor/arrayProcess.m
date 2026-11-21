function out = arrayProcess(arrayFiles, truthFile, dataDirOut, name, mode, opt)
%ARRAYPROCESS Anti-jamming array processing of the 7+1 CRPA element files.
%   out = arrayProcess(arrayFiles, truthFile, dataDirOut, name, mode, opt)
%
%   Reads the eight element files of generateIFdataArray.m and writes a data
%   set that the deep-tracking receiver reads unchanged:
%     <dataDirOut>\IFdata\navic_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin    (shared output)
%     <dataDirOut>\IFdata\navic_0IF_Fs5MHz_8bit_IQ_Fighter_<name>_prnNN.bin (per-satellite outputs)
%     <dataDirOut>\IFdata\..._time.mat, <dataDirOut>\FighterProfiles\ (truth, IMU)
%     <dataDirOut>\array_log_<name>.mat (beam phase, gain and SINR every 1 ms)
%   The receiver uses a per-satellite file for a NavIC channel when it exists.
%
%   mode
%     'single'  element 0 only (single-antenna reference)
%     'pi'      power inversion, one shared beam (conventional array unit)
%     'pi_ff'   'pi' plus per-satellite outputs whose carrier phase is
%               compensated for the beam response toward the INS-predicted
%               satellite direction, psi_i = arg(w^H a_i) (feed-forward, RQ4)
%     'mvdr'    per-satellite MVDR beams steered to the INS-predicted
%               directions (distortionless toward the model direction)
%   opt (optional)
%     tbMs        covariance/weight block [ms]: scalar, or 'sched' for the
%                 body-rate schedule T_b = min(Tmax, max(Tmin, kappa*dTheta/|w|))
%                 (default 'sched' for 'mvdr', 1 ms otherwise)
%     tbMinMs, tbMaxMs, kappa, nullWidthDeg   schedule (0.2, 10, 0.1, 20)
%     loadDb      diagonal loading relative to the noise power (default 0 dB)
%     attErrDeg   INS attitude error [roll pitch yaw] used for steering
%                 (default [0.01 -0.01 0.02], the receiver's initial error)
%     steerAtt    'ins' (default) or 'gnss': steer without INS attitude, with
%                 roll 0, pitch = flight-path angle, heading = ground track
%                 (no-INS ablation; combine with a fixed tbMs)
%     imuFile     IMU .mat (gyro for the schedule); default from the profile folder
%     profileDir  source FighterProfiles folder
%     outRmsCounts  output AGC target rms in 8-bit counts (default 25)
%     blankK      pulse blanking ahead of the spatial processing: samples whose
%                 amplitude on any element exceeds blankK times the noise rms
%                 are set to zero on all elements, with a guard of
%                 blankGuard samples (default: off; blankGuard 5)
%   Weights computed from block b are applied to block b+1 (causal).

if nargin < 6, opt = struct(); end
g = @(f, d) getfDefault(opt, f, d);
T = load(truthFile);
A = T.array;  M = A.M;  fs = 5e6;  nMs = 5000;
prn = T.prn(:)';  nSv = numel(prn);
nTotMs = floor(dir(arrayFiles{1}).bytes/(4*nMs));
isMvdr = strcmp(mode, 'mvdr');
tbOpt = g('tbMs', ternary(isMvdr, 'sched', 1));
tbMin = g('tbMinMs', 0.2);  tbMax = g('tbMaxMs', 10);  kappa = g('kappa', 0.1);  dTh = g('nullWidthDeg', 20);
loadLin = 10^(g('loadDb', 0)/10);
attErr = g('attErrDeg', [0.01 -0.01 0.02]);
steerAtt = g('steerAtt', 'ins');   % 'ins' (default) or 'gnss': no INS attitude, see below
profDir = g('profileDir', fileparts(fileparts(truthFile)));
targetRms = g('outRmsCounts', 25);
blankK = g('blankK', []);  blankGuard = g('blankGuard', 5);
noiseRms = [];  blankedFrac = 0;

% data set folders, timing file, truth/IMU
ifo = fullfile(dataDirOut, 'IFdata');  if ~exist(ifo, 'dir'), mkdir(ifo); end
pfo = fullfile(dataDirOut, 'FighterProfiles');  if ~exist(pfo, 'dir'), mkdir(pfo); end
base = ['navic_0IF_Fs5MHz_8bit_IQ_Fighter_' name];
[d0, b0] = fileparts(arrayFiles{1});
copyfile(fullfile(d0, [b0 '_time.mat']), fullfile(ifo, [base '_time.mat']));
for f = {['Fighter_' name '.csv'], ['Fighter_' name '_IMU_200Hz.mat']}
    src = fullfile(profDir, f{1});
    if exist(src, 'file'), copyfile(src, fullfile(pfo, f{1})); end
end
imuFile = g('imuFile', fullfile(profDir, ['Fighter_' name '_IMU_200Hz.mat']));
gyroRate = [];
if exist(imuFile, 'file')
    S = load(imuFile, 'imu');  gyroRate = struct('t', S.imu.t, 'w', sqrt(sum(S.imu.w.^2, 2)));   % rad/s
end

% outputs
fin = zeros(1, M);
for m = 1:M, fin(m) = fopen(arrayFiles{m}, 'r'); end
perSat = any(strcmp(mode, {'pi_ff', 'mvdr'}));
fShared = fopen(fullfile(ifo, [base '.bin']), 'w');
fSat = zeros(1, nSv);
if perSat
    for s = 1:nSv, fSat(s) = fopen(fullfile(ifo, sprintf('%s_prn%02d.bin', base, prn(s))), 'w'); end
end
cl = onCleanup(@() fcloseAll([fin fShared fSat]));

% steering directions: truth navigation-frame LOS (every 10 ms) rotated to the
% body frame with the attitude plus the INS attitude error
tRec = T.t(:);  attRec = T.att;  uN = T.uSatN;
Cerr = cnbFromEuler(attErr(1), attErr(2), attErr(3));
agcOut = ones(1, 1 + nSv)*NaN;
logEvery = 1;  nLog = floor(nTotMs/logEvery);   % 1 ms: the receiver feed-forward (utArrayAid) reads psi from this log
lg.t = zeros(nLog,1);  lg.tbMs = zeros(nLog,1);  lg.gainSat = zeros(nLog, nSv);  lg.psi = zeros(nLog, nSv);
lg.sinrGainDb = zeros(nLog, nSv);  lg.prn = prn;  lg.mode = mode;
lg.blankFrac = zeros(nLog, 1);  bfMs = 0;      % fraction of the samples of each ms blanked (receiver C/N0 option 'n')
w = repmat([1; zeros(M-1,1)], 1, nSv + 1);      % columns: shared, then per satellite
Racc = zeros(M);  accN = 0;  Rlast = eye(M);  sigma2 = 1;
% per-ms attitude and body rate, interpolated once (vectorised)
tMs = ((1:nTotMs)' - 0.5)*1e-3;
logTrue = strcmp(steerAtt, 'gnss');
if logTrue
    % ablation without INS attitude: roll 0, pitch = flight-path angle and
    % heading = ground track from the (error-free) GNSS velocity; the beam
    % response is also logged toward the true directions
    attTrueMs = interp1(tRec, attRec, min(max(((1:nTotMs)' - 0.5)*1e-3, tRec(1)), tRec(end)));
    attRec = gnssAttitude(tRec, T.posEcef);  Cerr = eye(3);
    lg.gainTrueDb = zeros(nLog, nSv);  lg.psiTrue = zeros(nLog, nSv);  lg.sinrGainTrueDb = zeros(nLog, nSv);
end
attMs = interp1(tRec, attRec, min(max(tMs, tRec(1)), tRec(end)));
if ~isempty(gyroRate), wMs = interp1(gyroRate.t, gyroRate.w, tMs, 'linear', 0); else, wMs = zeros(nTotMs, 1); end
% file I/O in batches of nBatch ms (one seek per file and batch, not per ms)
nBatch = 200;  readEl = 1:M;  if strcmp(mode, 'single'), readEl = 1; end
outBuf = cell(1, 1 + nSv);  for o = 1:1 + nSv, outBuf{o} = zeros(2*nMs*nBatch, 1, 'int8'); end
tStart = tic;
for k = 1:nTotMs
    kb = mod(k - 1, nBatch);                      % position within the batch
    if kb == 0
        nb = min(nBatch, nTotMs - k + 1);
        XB = zeros(M, nMs*nb);
        for m = readEl
            v = fread(fin(m), 2*nMs*nb, 'int16=>double');
            XB(m, :) = v(1:2:end) + 1i*v(2:2:end);
        end
    end
    tk = tMs(k);
    if strcmp(mode, 'single') && isempty(blankK)    % fast path: element 0 passed through (no blanking)
        Y = XB(1, kb*nMs + (1:nMs));
        p = mean(abs(Y).^2)/2;
        if isnan(agcOut(1)), agcOut(1) = targetRms/sqrt(p); else, agcOut(1) = 0.99*agcOut(1) + 0.01*targetRms/sqrt(p); end
        z = Y*agcOut(1);
        re = max(min(round(real(z)), 127), -127);  im = max(min(round(imag(z)), 127), -127);
        outBuf{1}(kb*2*nMs + (1:2*nMs)) = int8(reshape([re; im], [], 1));
        if kb == nb - 1, fwrite(fShared, outBuf{1}(1:2*nMs*nb), 'int8'); end
        continue
    end
    X = XB(:, kb*nMs + (1:nMs));
    % pulse blanking, identical on all elements (preserves the spatial structure)
    if ~isempty(blankK)
        % statistics over the elements read (single mode reads element 0 only;
        % the unread rows of X are zero and would set the noise estimate to 0)
        Xr = X(readEl, :);
        mag = max(abs(Xr), [], 1);
        % noise floor of this ms from the median magnitude before blanking
        % (Rayleigh median -> rms per component): robust while pulses occupy
        % less than half of the samples (DME: ~12 %), and it follows the onset of
        % a continuous jammer at once. A tracker that keeps the old floor blanks
        % everything after such an onset and never recovers.
        noiseRms = median(abs(Xr(:)))/sqrt(log(4));
        hit = mag > blankK*noiseRms;  bfMs = 0;
        if any(hit)
            hit = conv(double(hit), ones(1, 2*blankGuard + 1), 'same') > 0;
            X(:, hit) = 0;
            bfMs = mean(hit);  blankedFrac = blankedFrac + bfMs/nTotMs;
        end
    end
    % predicted body-frame LOS of each satellite at this ms
    att = attMs(k, :);
    Cnb = Cerr*cnbFromEuler(att(1), att(2), att(3));
    r = min(max(round(tk*100), 1), size(uN, 3));
    Asv = A.steer(Cnb*squeeze(uN(:, :, r)));
    if logTrue                                    % no-INS ablation: response toward the true directions
        attT = attTrueMs(k, :);
        AsvT = A.steer(cnbFromEuler(attT(1), attT(2), attT(3))*squeeze(uN(:, :, r)));
    end
    % block length for this ms
    if ischar(tbOpt)
        wr = wMs(k);
        tb = min(tbMax, max(tbMin, kappa*dTh*pi/180/max(wr, eps)*1e3));
        tb = min(tbMax, max(tbMin, 2^round(log2(tb))));                 % power-of-two grid
    else
        tb = tbOpt;
    end
    nbTarget = max(1, round(tb*nMs));                % samples per covariance block (tb in ms)
    Y = zeros(1 + nSv, nMs);
    pos = 1;
    while pos <= nMs
        % at least one sample: when the scheduled block shrinks (body rate up) the
        % samples already accumulated can exceed the new target; the block is then
        % closed after this sample
        len = min(max(nbTarget - accN, 1), nMs - pos + 1);
        idx = pos:pos + len - 1;
        Xb = X(:, idx);
        % apply the weights of the previous block (causal)
        switch mode
            case 'single', Y(1, idx) = Xb(1, :);
            case {'pi', 'pi_ff'}
                Y(1, idx) = w(:, 1)'*Xb;
                if perSat
                    for s = 1:nSv
                        psi = angle(w(:, 1)'*Asv(:, s));
                        Y(1 + s, idx) = Y(1, idx)*exp(-1i*psi);
                    end
                end
            case 'mvdr'
                Y(1, idx) = Xb(1, :);
                for s = 1:nSv, Y(1 + s, idx) = w(:, 1 + s)'*Xb; end
        end
        % accumulate the covariance; new weights when the block is complete
        Racc = Racc + Xb*Xb';  accN = accN + len;  pos = pos + len;
        if accN >= nbTarget && ~strcmp(mode, 'single')
            R = Racc/accN;  R = (R + R')/2;  Rlast = R;
            sigma2 = min(real(eig(R)));                % noise floor (fewer jammers than elements)
            Ri = inv(R + loadLin*sigma2*eye(M));
            switch mode
                case {'pi', 'pi_ff'}
                    w(:, 1) = Ri(:, 1)/Ri(1, 1);
                case 'mvdr'
                    for s = 1:nSv
                        q = Ri*Asv(:, s);  w(:, 1 + s) = q/(Asv(:, s)'*q);
                    end
            end
            Racc = zeros(M);  accN = 0;
        elseif accN >= nbTarget
            Racc = zeros(M);  accN = 0;
        end
    end
    % output AGC (slow, per stream) and 8-bit quantization
    for o = 1:(1 + nSv*perSat)
        p = mean(abs(Y(o, :)).^2)/2;
        if isnan(agcOut(o)), agcOut(o) = targetRms/sqrt(p); else, agcOut(o) = 0.99*agcOut(o) + 0.01*targetRms/sqrt(p); end
        z = Y(o, :)*agcOut(o);
        re = max(min(round(real(z)), 127), -127);  im = max(min(round(imag(z)), 127), -127);
        outBuf{o}(kb*2*nMs + (1:2*nMs)) = int8(reshape([re; im], [], 1));
        if kb == nb - 1                           % end of the batch: write it
            if o == 1, fo = fShared; else, fo = fSat(o - 1); end
            fwrite(fo, outBuf{o}(1:2*nMs*nb), 'int8');
        end
    end
    % diagnostics
    if mod(k, logEvery) == 0
        i = k/logEvery;  lg.t(i) = k*1e-3;  lg.tbMs(i) = tb;  lg.blankFrac(i) = bfMs;
        for s = 1:nSv
            if strcmp(mode, 'mvdr'), ws = w(:, 1 + s); else, ws = w(:, 1); end
            resp = ws'*Asv(:, s);
            lg.gainSat(i, s) = 10*log10(abs(resp)^2/real(ws'*ws));
            lg.psi(i, s) = angle(resp);
            lg.sinrGainDb(i, s) = 10*log10(abs(resp)^2*sigma2/real(ws'*Rlast*ws));   % output SINR relative to one interference-free element
            if logTrue
                rT = ws'*AsvT(:, s);
                lg.gainTrueDb(i, s) = 10*log10(abs(rT)^2/real(ws'*ws));
                lg.psiTrue(i, s) = angle(rT);
                lg.sinrGainTrueDb(i, s) = 10*log10(abs(rT)^2*sigma2/real(ws'*Rlast*ws));
            end
        end
    end
    if mod(k, 10000) == 0, fprintf('  %s: %d ms of %d (%.1f min)\n', mode, k, nTotMs, toc(tStart)/60); end
end
lg.blankedFraction = blankedFrac;
save(fullfile(dataDirOut, ['array_log_' name '.mat']), '-struct', 'lg');
out = struct('dataDir', dataDirOut, 'minutes', toc(tStart)/60, 'mode', mode, 'blankedFraction', blankedFrac);
fprintf('Array processing (%s) done in %.1f min -> %s\n', mode, out.minutes, dataDirOut);
end

function att = gnssAttitude(t, pe)
% [roll pitch heading] in degrees from the ECEF track: roll 0, pitch = flight-path
% angle, heading = course over ground (velocity by central differences)
v = [gradient(pe(:,1), t), gradient(pe(:,2), t), gradient(pe(:,3), t)];
a = 6378137;  e2 = 6.69437999014e-3;
lon = atan2(pe(:,2), pe(:,1));  p = hypot(pe(:,1), pe(:,2));
lat = atan2(pe(:,3), p*(1 - e2));
for it = 1:5
    N = a./sqrt(1 - e2*sin(lat).^2);
    h = p./cos(lat) - N;
    lat = atan2(pe(:,3), p.*(1 - e2*N./(N + h)));
end
sl = sin(lat); cl = cos(lat); so = sin(lon); co = cos(lon);
vn = -sl.*co.*v(:,1) - sl.*so.*v(:,2) + cl.*v(:,3);
ve = -so.*v(:,1) + co.*v(:,2);
vd = -cl.*co.*v(:,1) - cl.*so.*v(:,2) - sl.*v(:,3);
att = [zeros(size(t)), atan2d(-vd, hypot(vn, ve)), mod(atan2d(ve, vn), 360)];
end

function fcloseAll(f)
for x = f(f > 0), fclose(x); end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end

function v = getfDefault(S, f, d)
if isfield(S, f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end

function C = cnbFromEuler(rollD, pitchD, headD)
sp = sind(rollD); cp = cosd(rollD); st = sind(pitchD); ct = cosd(pitchD); ss = sind(headD); cs = cosd(headD);
C = [ct*cs, ct*ss, -st; -cp*ss + sp*st*cs, cp*cs + sp*st*ss, sp*ct; sp*ss + cp*st*cs, -sp*cs + cp*st*ss, cp*ct];
end
