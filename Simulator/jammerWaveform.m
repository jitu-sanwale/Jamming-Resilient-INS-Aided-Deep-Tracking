function [x, st] = jammerWaveform(st, n, fs, t0)
%JAMMERWAVEFORM Next n complex baseband samples (unit power) of one jammer.
%   [x, st] = jammerWaveform(st, n, fs, t0)
%
%   st   - jammer state; on the first call a struct with field type and the
%          parameters below (missing ones take the defaults); returned updated
%   n    - number of samples;  fs - sample rate [Hz];  t0 - time of the first sample [s]
%
%   Types and parameters (offsets from the NavIC L5 carrier):
%     'bbn'    band-limited Gaussian noise: bandwidthHz (default 2.046e6, the
%              BPSK(1) main lobe), offsetHz (0). Unit average power.
%     'cw'     tone: offsetHz (default 0). Unit power.
%     'chirp'  linear sweep: f0Hz (-2e6), f1Hz (+2e6), periodS (0.5e-3), as the
%              PPD-type model of ChirpJammerSetting.m. Unit power.
%     'dme'    one DME/TACAN beacon: Gaussian pulse pairs exp(-a t^2/2) +
%              exp(-a (t - dt)^2/2), a = 4.5e11 s^-2 (half-amplitude width
%              3.5 us), dt = spacingS (12e-6 X mode; 30e-6/36e-6 Y mode), Poisson
%              pulse-pair arrivals at rateHz (DME 700-2700 ppps, TACAN up to
%              3600; default 2700), carrier offset offsetHz (+0.55e6, i.e. the
%              1177-MHz channel; channels lie on a 1-MHz grid). Parameters from
%              Lestarquit et al. (NAVITEC 2022), Garcia-Pena et al. (ITM 2020),
%              Musumeci et al. (IET RSN 2014). Unit PEAK power; a hotspot is
%              modelled as several beacons, each a separate jammer.
%   Phase continuity across calls is kept in st (filter state, phase, pulse times).

if ~isfield(st, 'init') || ~st.init
    st = initState(st, fs);
end
t = t0 + (0:n-1)/fs;
switch st.type
    case 'bbn'
        w = (randn(1, n) + 1i*randn(1, n))/sqrt(2);
        [y, st.zf] = filter(st.b, 1, w, st.zf);
        x = y/st.gainNorm .* exp(1i*2*pi*st.offsetHz*t);
    case 'cw'
        x = exp(1i*(2*pi*st.offsetHz*t + st.phi0));
    case 'chirp'
        tl = mod(t, st.periodS);
        k = (st.f1Hz - st.f0Hz)/st.periodS;
        x = exp(1i*(2*pi*(st.f0Hz*tl + 0.5*k*tl.^2) + st.phi0));
    case 'dme'
        % draw pulse-pair arrival times up to the end of this block
        tEnd = t(end) + 5*st.widthS + st.spacingS;
        while st.nextT < tEnd
            st.pairs(end+1) = st.nextT;                      %#ok<AGROW>
            st.nextT = st.nextT - log(rand)/st.rateHz;       % Poisson arrivals
        end
        st.pairs(st.pairs < t(1) - 5*st.widthS - st.spacingS) = [];
        sig = st.widthS/(2*sqrt(2*log(2)));                 % half-amplitude width: exp(-(w/2)^2/(2 sig^2)) = 0.5
        env = zeros(1, n);
        for p = st.pairs
            env = env + exp(-(t - p).^2/(2*sig^2)) + exp(-(t - p - st.spacingS).^2/(2*sig^2));
        end
        x = env .* exp(1i*(2*pi*st.offsetHz*t + st.phi0));
    otherwise
        error('jammerWaveform: unknown type %s', st.type);
end
end

function st = initState(st, fs)
d = struct('bandwidthHz', 2.046e6, 'offsetHz', 0, 'f0Hz', -2e6, 'f1Hz', 2e6, 'periodS', 0.5e-3, ...
           'rateHz', 2700, 'spacingS', 12e-6, 'widthS', 3.5e-6);
if strcmp(st.type, 'dme') && ~isfield(st, 'offsetHz'), d.offsetHz = 0.55e6; end
f = fieldnames(d);
for k = 1:numel(f)
    if ~isfield(st, f{k}) || isempty(st.(f{k})), st.(f{k}) = d.(f{k}); end
end
st.phi0 = 2*pi*rand;
if strcmp(st.type, 'bbn')
    bw = min(st.bandwidthHz, 0.98*fs);
    st.b = fir1(64, bw/fs);                 % complex baseband: two-sided bandwidth bw
    st.zf = zeros(numel(st.b) - 1, 1);
    st.gainNorm = sqrt(sum(abs(st.b).^2));  % unit output power for unit-power white input
end
if strcmp(st.type, 'dme')
    st.pairs = [];
    st.nextT = -log(rand)/st.rateHz;
end
st.init = true;
end
