function out = generateIFdataArray(name, tag, opt)
%GENERATEIFDATAARRAY NavIC L5 IF files of the 7+1 CRPA for one trajectory and scenario.
%   out = generateIFdataArray(name, tag, opt)
%
%   name - profile name; trajectory <opt.profileDir>\Fighter_<name>.csv
%   tag  - scenario tag used in the file names (e.g. 'clean', 'S1_seed03')
%   opt  - struct:
%     profileDir   trajectory folder (default <simulator>\FighterProfiles)
%     outDir       output folder (default <simulator>\IFarray)
%     ms           record length [ms] (default 60000)
%     seed         noise seed (default 0)
%     artificial   true: seven-satellite constellation (default true)
%     jammers      struct array (default none), fields:
%                    type   'bbn' | 'cw' | 'chirp' | 'dme' (jammerWaveform.m)
%                    pos    [lat lon h] (deg, deg, m) or @(t) -> [lat lon h]
%                    jsDb   J/S [dB] at an isotropic element w.r.t. a 45 dB-Hz
%                           NavIC signal: scalar or [t js] ramp table
%                    erpDbw (optional, instead of jsDb) transmitter EIRP [dBW]
%                           (peak EIRP for DME, 37 dBW typical); the J/S then
%                           follows from free-space loss and N0 = -201.5 dBW/Hz
%                    tOn, tOff  active window [s from the start]
%                    params struct of waveform parameters (may be empty)
%     adcBits      ADC resolution (default 14)
%     agcBackoffDb AGC target rms below full scale (default 12 dB)
%     agcTauMs     AGC time constant [ms] (default 10)
%     couplingDb   adjacent-element coupling (default -20; -Inf: none)
%     cnoDb        satellite C/N0 [dB-Hz] before antenna gains (default 45, as initSettings)
%   out - struct: files (8 element files), truthFile, info.
%
%   Files: <outDir>\navic8_e<m>_0IF_Fs5MHz_16bit_IQ_Fighter_<name>_<tag>.bin,
%   m = 0..7 (0: centre/reference element), int16 I/Q at 5 MHz, and
%   <outDir>\navic8_truth_Fighter_<name>_<tag>.mat (geometry, attitude, AGC).

simDir = fileparts(mfilename('fullpath'));
if nargin < 3, opt = struct(); end
g = @(f, d) getfDefault(opt, f, d);
profDir = g('profileDir', fullfile(simDir, 'FighterProfiles'));
outDir  = g('outDir', fullfile(simDir, 'IFarray'));
if ~exist(outDir, 'dir'), mkdir(outDir); end
trajFile = fullfile(profDir, ['Fighter_' name '.csv']);
if ~exist(trajFile, 'file'), error('generateIFdataArray: trajectory %s not found.', trajFile); end

here = pwd;  back = onCleanup(@() cd(here));
cd(simDir);
settings = initSettings();
settings.trajFile = trajFile;
settings.gnss = 'NavIC';
settings.msToProcess = g('ms', 60000);
settings.defaultCNoDB = g('cnoDb', settings.defaultCNoDB);   % satellite C/N0 [dB-Hz] (default 45; tests use more)
if g('artificial', true), settings.navicArtificial = navicArtificialDefault(); end

cfg.array = arrayGeometry(struct('couplingDb', g('couplingDb', -20)));
cfg.jammers = g('jammers', struct('type', {}, 'pos', {}, 'jsDb', {}, 'tOn', {}, 'tOff', {}, 'params', {}));
cfg.adcBits = g('adcBits', 14);
cfg.agcBackoffDb = g('agcBackoffDb', 12);
cfg.agcAlpha = 1/g('agcTauMs', 10);              % AGC smoothing per 1-ms block
cfg.files = cell(1, cfg.array.M);
for m = 1:cfg.array.M
    cfg.files{m} = fullfile(outDir, sprintf('navic8_e%d_0IF_Fs5MHz_16bit_IQ_Fighter_%s_%s.bin', m - 1, name, tag));
end
cfg.truthFile = fullfile(outDir, sprintf('navic8_truth_Fighter_%s_%s.mat', name, tag));

fprintf('=== Array IF for %s / %s, %.1f s, %d jammer(s) -> %s\n', name, tag, settings.msToProcess/1000, ...
    numel(cfg.jammers), outDir);
rngState = rng;  rng(g('seed', 0), 'twister');
info = signalGeneratingArray(settings, cfg);
rng(rngState);
out = struct('files', {cfg.files}, 'truthFile', cfg.truthFile, 'info', info);
end

function v = getfDefault(S, f, d)
if isfield(S, f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
