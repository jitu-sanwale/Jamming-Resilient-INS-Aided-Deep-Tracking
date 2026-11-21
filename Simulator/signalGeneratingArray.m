function info = signalGeneratingArray(settings, cfg)
%SIGNALGENERATINGARRAY NavIC L5 IF samples of the 7+1 CRPA, with jammers.
%   info = signalGeneratingArray(settings, cfg)
%
%   Multi-element version of signalGenerating.m (NavIC L5 SPS path only; no
%   multipath or front-end filter). Per 1-ms block and satellite the signal is
%   generated once at the array centre, as in the single-antenna simulator
%   (travel time, satellite clock, IRNWT time model, receiver L5 delay,
%   ionosphere and troposphere, Rx and Tx antenna gain, NavIC navigation
%   message), and then given to each element m the carrier phase of its
%   displacement d_m towards the satellite:
%       x_m = g * [C * exp(-1i*w_c*(u'*d)/c)]_m * s,
%   with u the receiver-to-satellite unit vector and C the coupling matrix
%   (arrayGeometry.m). The code delay across the 25-cm aperture (< 0.9 ns, below
%   5e-3 of a chip) is neglected (narrowband array model, as the paper).
%   Element gains: the attitude-dependent Rx pattern of the single-antenna
%   simulator, applied to every (identical, zenith-pointing) element.
%   Jammers (cfg.jammers, see generateIFdataArray.m) are added per element from
%   their geometric direction and range, with a carrier phase that follows the
%   range (Doppler). Receiver noise is independent per element.
%   One AGC gain, common to all elements, follows the total input power with a
%   time constant of 1/cfg.agcAlpha ms (default 10 ms); samples are quantized to cfg.adcBits bits (clipping)
%   and written as int16 I/Q, one file per element.
%
%   cfg fields: files (1x8 cellstr, output IF files), truthFile, array (from
%   arrayGeometry), jammers (struct array, may be empty), adcBits, agcBackoffDb.
%   info: struct with the PRN list, timing and AGC statistics.

A = cfg.array;  M = A.M;
blockSize = round(settings.samplingFreq*0.001);
blockTime = blockSize/settings.samplingFreq;
c = settings.c;

%% trajectory (as signalGenerating.m)
tr = readmatrix(settings.trajFile, 'NumHeaderLines', 1);
if settings.msToProcess/1000 > tr(end,1), settings.msToProcess = floor(tr(end,1)*1000); end
iterCnt = round(settings.msToProcess/1000/blockTime);
qTime = min((0:iterCnt)*blockTime, tr(end,1));
q = interp1(tr(:,1), tr(:,2:7), qTime, 'linear');            % lat lon alt roll pitch heading
traj = zeros(iterCnt + 1, 9);
for k = 1:iterCnt + 1
    [x, y, z] = geo2cartd(q(k,1), q(k,2), q(k,3), 5);
    traj(k,:) = [q(k,1:3), x, y, z, q(k,4:6)];
end

%% ephemerides (GPS time base, NavIC records, optional artificial satellites)
[ephG, ionoutc] = rinexeV3(settings.rinexfile);
ionoAlpha = [ionoutc.alpha0, ionoutc.alpha1, ionoutc.alpha2, ionoutc.alpha3];
ionoBeta  = [ionoutc.beta0,  ionoutc.beta1,  ionoutc.beta2,  ionoutc.beta3];
startTime = ephG(1).toc;  gpsWeek = ephG(1).weekNrm;
eph = rinexeIRNSS(settings.rinexfile, startTime);
if isfield(settings, 'navicArtificial') && ~isempty(settings.navicArtificial)
    eph = navicArtificialEph(eph, startTime, settings.navicArtificial);
end
settings.carrFreqBasis = A.fc;
igt = struct('A0', 25e-9, 'A1', 2e-14, 'A2', 0, 'tot', 16*floor((startTime - 3600)/16), 'WNot', gpsWeek - 1024);
navicDtErr = 3e-9;  navicRxL5Delay = 12e-9;
navicDtIG = @(tg) igt.A0 + igt.A1*check_t(tg - igt.tot) + igt.A2*check_t(tg - igt.tot).^2 + navicDtErr;
navicTime = struct('igt', igt, 'dtErr', navicDtErr, 'rxL5Delay', navicRxL5Delay, ...
    'isbTruth_m', c*(navicRxL5Delay - navicDtErr), 'startTime', startTime); %#ok<NASGU>
for m = 1:M
    [d, b] = fileparts(cfg.files{m});
    save(fullfile(d, [b '_time.mat']), '-struct', 'navicTime');
end
ionoScale = (1575.42e6/A.fc)^2;

%% visible satellites, amplitudes, patterns, messages, codes
[satList, el0, az0] = getVisibleSat(eph, startTime, traj(1,4:6)', settings);
nSv = numel(satList);
fprintf('Array IF: %d NavIC satellites visible: %s\n', nSv, mat2str(satList));
CNo = 10.^((CNoSetting(satList, settings.defaultCNoDB) + 1.5)/10);
carrAmp = 2*settings.gwnAmp*sqrt(CNo/settings.samplingFreq)/sqrt(2);      % I/Q file (fileType 2)
cRef = (2*settings.gwnAmp*sqrt(10^4.5/settings.samplingFreq)/sqrt(2))^2;  % carrier power of a 45 dB-Hz signal
rxPat = loadAntennaPattern(settings.rxAntennaFile);
txPat = loadAntennaPattern(settings.txAntennaFile);
for s = 1:nSv
    prn = satList(s);
    navBits(s).template = navic_eph2sbf(eph(prn), ionoutc, igt); %#ok<AGROW>
    navBits(s).refTime = -Inf;  navBits(s).dataBit = [];
    code = generateCAcode_L5_band(prn);
    codeTab(s,:) = [code(end) code code(1)]; %#ok<AGROW>
end
bitPeriod = 0.02;
localOs = (A.fc - settings.IF)*2*pi;  carrRad = A.fc*2*pi;

%% jammers
jam = cfg.jammers;  nJ = numel(jam);
for j = 1:nJ
    if ~isfield(jam(j), 'state') || isempty(jam(j).state), jam(j).state = struct('type', jam(j).type); end
    pr = jam(j).params;  fn = fieldnames(pr);
    for f = 1:numel(fn), jam(j).state.(fn{f}) = pr.(fn{f}); end
end

%% output files and truth record (every 10 ms)
fid = zeros(1, M);
for m = 1:M
    fid(m) = fopen(cfg.files{m}, 'w');
    if fid(m) < 0, error('signalGeneratingArray: cannot open %s', cfg.files{m}); end
end
cleanup = onCleanup(@() arrayfun(@(f) fclose(f), fid(fid > 0)));
nRec = floor(iterCnt/10);
truth.t = zeros(nRec,1);  truth.att = zeros(nRec,3);  truth.posEcef = zeros(nRec,3);
truth.prn = satList;  truth.uSatB = zeros(3, nSv, nRec);  truth.uSatN = zeros(3, nSv, nRec);
truth.uJamB = zeros(3, nJ, nRec);  truth.jsDb = zeros(nJ, nRec);  truth.agcGain = zeros(nRec,1);
truth.array = A;  truth.posB = A.posB;  truth.C = A.C;
fullScale = 2^(cfg.adcBits - 1) - 1;
targetRms = fullScale/10^(cfg.agcBackoffDb/20);
agc = [];  clipCount = 0;

%% main loop
RxTime = [startTime, startTime + blockTime];
RxPos = traj(1:2, 4:6);
fprintf('Array IF: generating %d ms, %d elements ...\n', iterCnt, M);
tGen = tic;
for k = 1:iterCnt
    llh = traj(k,1:3);  rollD = traj(k,7);  pitchD = traj(k,8);  headD = traj(k,9);
    Cnb = cnbFromEuler(rollD, pitchD, headD);           % NED -> body
    Cen = cenMat(llh(1), llh(2));                       % ECEF -> NED
    Ceb = Cnb*Cen;                                      % ECEF -> body
    X = complex(zeros(M, blockSize));
    rxS = linspace(RxTime(1), RxTime(2) - blockTime/blockSize, blockSize);
    rec = mod(k, 10) == 0;  r = k/10;

    for s = 1:nSv
        prn = satList(s);
        [Tx, svClk, svPos] = GetTravelTime(RxTime, RxPos, eph(prn), settings);
        Tx = Tx + svClk;
        Tx = Tx + navicDtIG(Tx) - navicRxL5Delay;
        [azRx, elRx] = topocent(RxPos(1,:)', (svPos(1,:) - RxPos(1,:))');
        [azB, elB] = applyAttitudeToBodyFrame(azRx, elRx, rollD, pitchD, headD);
        gDb = lookupAntennaGain(rxPat, azB, elB, settings.rxBelowHorizonSlopeDb, settings.rxBelowHorizonFloorDb);
        toC = -svPos(1,:);  toR = RxPos(1,:) - svPos(1,:);
        gDb = gDb + lookupAntennaGain(txPat, 0, acosd(dot(toC, toR)/(norm(toC)*norm(toR))));
        atm = ionoScale*klobuchar(ionoAlpha, ionoBeta, RxTime(1), llh(1), llh(2), elRx, azRx) ...
              + saastamoinenTropo(elRx, llh(1), llh(3), settings.tropoP0, settings.tropoT0, settings.tropoRH);
        Tx = Tx - atm;
        txS = linspace(Tx(1), Tx(2) - (Tx(2) - Tx(1))/blockSize, blockSize);
        carr = exp(-1i*(carrRad*txS - localOs*rxS));
        code = codeTab(s, ceil(mod(txS*settings.codeFreqBasis, settings.codeLength)) + 1);
        % NavIC symbols (as signalGenerating.m)
        while Tx(end) >= navBits(s).refTime + 48
            if navBits(s).refTime == -Inf, seedT = Tx(1); else, seedT = navBits(s).refTime + 48 + 0.5; end
            [sym, newRef] = navic_generateNavMsg(navBits(s).template, seedT);
            navBits(s).prevBit = navBits(s).dataBit;  navBits(s).prevRef = navBits(s).refTime;
            navBits(s).dataBit = reshape(sym.', 1, 2400);  navBits(s).refTime = newRef;
        end
        symIdx = floor((txS - navBits(s).refTime)/bitPeriod) + 1;
        bit = zeros(1, blockSize);  inC = symIdx >= 1;
        bit(inC) = navBits(s).dataBit(symIdx(inC));
        if any(~inC), bit(~inC) = navBits(s).prevBit(floor((txS(~inC) - navBits(s).prevRef)/bitPeriod) + 1); end
        sig = (code .* bit .* carr) * carrAmp(s) * 10^(gDb/20);
        % element phases from the array geometry (body frame)
        uE = (svPos(1,:) - RxPos(1,:))';  uE = uE/norm(uE);
        uB = Ceb*uE;
        X = X + A.steer(uB) * sig;
        if rec, truth.uSatB(:,s,r) = uB;  truth.uSatN(:,s,r) = Cen*uE; end
    end

    % jammers
    tEl = RxTime(1) - startTime;
    for j = 1:nJ
        if tEl < jam(j).tOn || tEl >= jam(j).tOff, continue; end
        pJ0 = jammerPosEcef(jam(j), tEl);  pJ1 = jammerPosEcef(jam(j), tEl + blockTime);
        dv0 = pJ0 - RxPos(1,:)';  dv1 = pJ1 - RxPos(2,:)';
        R0 = norm(dv0);  R1 = norm(dv1);
        uE = dv0/R0;  uB = Ceb*uE;
        [azJ, elJ] = topocent(RxPos(1,:)', dv0);
        [azJB, elJB] = applyAttitudeToBodyFrame(azJ, elJ, rollD, pitchD, headD);
        gJ = lookupAntennaGain(rxPat, azJB, elJB, settings.rxBelowHorizonSlopeDb, settings.rxBelowHorizonFloorDb);
        jsDb = jammerJS(jam(j), tEl, R0, settings.c/A.fc);
        [w, jam(j).state] = jammerWaveform(jam(j).state, blockSize, settings.samplingFreq, tEl);
        Rs = R0 + (R1 - R0)*(0:blockSize-1)/blockSize;               % range across the block
        js = sqrt(cRef*10^(jsDb/10)) * 10^(gJ/20) * w .* exp(-1i*carrRad*Rs/c);
        X = X + A.steer(uB) * js;
        if rec, truth.uJamB(:,j,r) = uB;  truth.jsDb(j,r) = jsDb; end
    end

    % receiver noise, independent per element
    X = X + settings.gwnAmp*(randn(M, blockSize) + 1i*randn(M, blockSize));

    % common AGC (10-ms time constant on the total power, as a hardware AGC
    % that settles within milliseconds of a jammer onset) and ADC
    p = mean(abs(X(:)).^2)/2;                           % per-component power
    if isempty(agc), agc = targetRms/sqrt(p); else, agc = agc*(1 - cfg.agcAlpha) + cfg.agcAlpha*targetRms/sqrt(p); end
    Y = X*agc;
    re = round(real(Y));  im = round(imag(Y));
    clipCount = clipCount + sum(abs(re(:)) > fullScale) + sum(abs(im(:)) > fullScale);
    re = max(min(re, fullScale), -fullScale);  im = max(min(im, fullScale), -fullScale);
    for m = 1:M
        fwrite(fid(m), reshape([re(m,:); im(m,:)], [], 1), 'int16');
    end

    if rec
        truth.t(r) = tEl + blockTime;  truth.att(r,:) = [rollD pitchD headD];
        truth.posEcef(r,:) = RxPos(1,:);  truth.agcGain(r) = agc;
    end
    RxTime = RxTime + blockTime;
    RxPos(1,:) = RxPos(2,:);  RxPos(2,:) = traj(min(k + 2, iterCnt + 1), 4:6);
    if mod(k, 5000) == 0
        fprintf('  %6d ms of %d  (%.1f min, clipped %.2e of samples)\n', k, iterCnt, toc(tGen)/60, ...
            clipCount/(2*M*blockSize*k));
    end
end
truth.clipFraction = clipCount/(2*M*blockSize*iterCnt);
if isfield(jam, 'state'), jam = rmfield(jam, 'state'); end
truth.jammers = jam;
truth.settings = settings;  truth.az0 = az0;  truth.el0 = el0; %#ok<STRNU>
save(cfg.truthFile, '-struct', 'truth', '-v7.3');
info = struct('satList', satList, 'minutes', toc(tGen)/60, 'clipFraction', truth.clipFraction);
fprintf('Array IF done in %.1f min; clipped fraction %.2e\n', info.minutes, info.clipFraction);
end

%% ------------------------------------------------------------------------
function p = jammerPosEcef(j, t)
% jammer position (ECEF, column) at elapsed time t: fixed [lat lon h] or a function handle
if isa(j.pos, 'function_handle'), lla = j.pos(t); else, lla = j.pos; end
[x, y, z] = geo2cartd(lla(1), lla(2), lla(3), 5);
p = [x; y; z];
end

function v = jammerJS(j, t, R, lam)
% J/S [dB] at an isotropic reference element, relative to a 45 dB-Hz NavIC
% signal. Either given directly (jsDb: scalar or [t js] ramp table), or from
% the transmitter EIRP and free-space loss (erpDbw: scalar or [t erp] table),
% with the receiver noise density N0 = -201.5 dBW/Hz (Lestarquit et al., 2022;
% Garcia-Pena et al., 2020), i.e. C(45 dB-Hz) = -156.5 dBW. For DME, erpDbw
% is the PEAK EIRP (37 dBW typical), so v is the peak J/S.
if isfield(j, 'erpDbw') && ~isempty(j.erpDbw)
    if isscalar(j.erpDbw), erp = j.erpDbw; else, erp = interp1(j.erpDbw(:,1), j.erpDbw(:,2), t, 'linear', 'extrap'); end
    v = erp - 20*log10(4*pi*R/lam) - (-201.5 + 45);
elseif isscalar(j.jsDb)
    v = j.jsDb;
else
    v = interp1(j.jsDb(:,1), j.jsDb(:,2), t, 'linear', 'extrap');
end
end

function C = cnbFromEuler(rollD, pitchD, headD)
% NED -> body, as applyAttitudeToBodyFrame.m
sp = sind(rollD); cp = cosd(rollD); st = sind(pitchD); ct = cosd(pitchD); ss = sind(headD); cs = cosd(headD);
C = [ct*cs, ct*ss, -st; -cp*ss + sp*st*cs, cp*cs + sp*st*ss, sp*ct; sp*ss + cp*st*cs, -sp*cs + cp*st*ss, cp*ct];
end

function C = cenMat(latD, lonD)
% ECEF -> NED
sL = sind(latD); cL = cosd(latD); sl = sind(lonD); cl = cosd(lonD);
C = [-sL*cl, -sL*sl, cL; -sl, cl, 0; -cL*cl, -cL*sl, -sL];
end
