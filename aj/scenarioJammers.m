function [jam, info] = scenarioJammers(scen, seed)
%SCENARIOJAMMERS Jammer sets of the navigation-level scenarios S1-S3 (paper Table VIII).
%   [jam, info] = scenarioJammers(scen, seed)
%   scen - 'S1' | 'S2' | 'S3';  seed - Monte Carlo index (geometry draw)
%   jam  - struct array for generateIFdataArray (opt.jammers)
%   info - struct describing the draw (for the result tables)
%
%   All profiles start at 20.1085 N, 73.9532 E, 5000 m, flying east at 250 m/s
%   (the baseline's first 12 s are common to all profiles). Ground jammers are
%   at h = 0, i.e. below the aircraft's horizon.
%   S1  broadband jamming: 1-3 BBN ground jammers, ground range 15-50 km, any
%       bearing; on at 15 s, J/S of each ramped from 40 dB (15 s) to 100 dB
%       (45 s), off at 50 s (tolerance, coasting, recovery).
%   S2  jammer near a NavIC satellite: one airborne BBN jammer flying with the
%       aircraft (250 m/s east) 8-12 deg below the direction of NavIC PRN 10
%       (az 101.2, el 24.8 deg), J/S 70 dB from 15 s to 50 s.
%   S3  compound L5 jamming: four DME beacons on the 1175-1178 MHz channels
%       (peak EIRP 37 dBW, 700-2700 ppps, X mode) at 30-150 km on the ground,
%       a PPD-type chirp jammer (J/S 50 dB) and one BBN jammer (J/S 60 dB) on
%       the ground at 10-30 km; chirp and BBN on from 15 s to 50 s, DME always.
rng(1000*double(scen(2)) + seed, 'twister');
lat0 = 20.1085;  lon0 = 73.9532;  h0 = 5000;
geo = @(rngKm, brgDeg, h) [lat0 + rngKm*cosd(brgDeg)/111.32, ...
                           lon0 + rngKm*sind(brgDeg)/(111.32*cosd(lat0)), h];
blank = struct('type', {}, 'pos', {}, 'jsDb', {}, 'erpDbw', {}, 'tOn', {}, 'tOff', {}, 'params', {});
jam = blank;  info = struct('scen', scen, 'seed', seed);
switch scen
    case 'S1'
        n = randi(3);
        for k = 1:n
            jam(k).type = 'bbn';  jam(k).pos = geo(15 + 35*rand, 360*rand, 0);
            jam(k).jsDb = [15 40; 45 100];  jam(k).erpDbw = [];
            jam(k).tOn = 15;  jam(k).tOff = 50;  jam(k).params = struct();
        end
        info.nJammers = n;
    case 'S2'
        sep = 8 + 4*rand;                          % deg below PRN 10
        az = 101.2;  el = 24.8 - sep;  dh = 7000;  % jammer at 12 km altitude
        dist = dh/tand(el)/1e3;                    % ground distance [km]
        p0 = geo(dist, az, h0 + dh);
        vEast = 250/(111.32e3*cosd(lat0));         % deg/s, flies with the aircraft
        jam(1).type = 'bbn';  jam(1).pos = @(t) [p0(1), p0(2) + vEast*t, p0(3)];
        jam(1).jsDb = 70;  jam(1).erpDbw = [];  jam(1).tOn = 15;  jam(1).tOff = 50;  jam(1).params = struct();
        info.sepDeg = sep;
    case 'S3'
        chans = [-1.45 -0.45 0.55 1.55]*1e6;      % 1175-1178 MHz relative to 1176.45 MHz
        rates = [2700 1350 2700 700];
        for k = 1:4
            jam(k).type = 'dme';  jam(k).pos = geo(30 + 120*rand, 360*rand, 0);
            jam(k).jsDb = [];  jam(k).erpDbw = 37;  jam(k).tOn = 0;  jam(k).tOff = Inf;
            jam(k).params = struct('offsetHz', chans(k), 'rateHz', rates(k));
        end
        jam(5).type = 'chirp';  jam(5).pos = geo(10 + 20*rand, 360*rand, 0);
        jam(5).jsDb = 50;  jam(5).erpDbw = [];  jam(5).tOn = 15;  jam(5).tOff = 50;  jam(5).params = struct();
        jam(6).type = 'bbn';  jam(6).pos = geo(10 + 20*rand, 360*rand, 0);
        jam(6).jsDb = 60;  jam(6).erpDbw = [];  jam(6).tOn = 15;  jam(6).tOff = 50;  jam(6).params = struct();
    otherwise
        error('scenarioJammers: unknown scenario %s', scen);
end
end
