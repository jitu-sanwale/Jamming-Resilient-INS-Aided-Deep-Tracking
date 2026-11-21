function eph = navicArtificialEph(eph, refSow, sats)
%NAVICARTIFICIALEPH Add artificial NavIC satellites as broadcast-ephemeris records.
%   eph = navicArtificialEph(eph, refSow, sats)
%
%   eph    - NavIC eph(prn) struct array from rinexeIRNSS.m (real satellites)
%   refSow - reference time, seconds of week (the simulation start time)
%   sats   - struct array, one element per artificial satellite, fields
%              prn          NavIC PRN (1-14, not used by a real satellite)
%              inclDeg      inclination [deg] (0 for GEO)
%              lonCrossDeg  longitude [deg E] of the equator crossing (GEO slot)
%              u0Deg        argument of latitude at refSow [deg]
%                           (90: northern lobe, 270: southern lobe; 0 for GEO)
%            default: navicArtificialDefault()
%
%   Each satellite gets a circular geosynchronous orbit in the broadcast
%   (IS-GPS-200 / IRNSS ICD) Keplerian form read by getSatPos.m and encoded by
%   navic_eph2sbf.m: sqrt(A) from a mean motion equal to the Earth rotation
%   rate, e = 0, omega = 0, deltan = Omegadot = idot = 0, no harmonic
%   corrections, zero clock terms and group delay, healthy. With n = OmegaE,
%   the Earth-fixed track is the GEO slot (i = 0) or the IGSO figure eight
%   crossing the equator at lonCrossDeg (u = 0 there):
%     lon(t) = lonCross - u(t) + atan2(cos i sin u(t), cos u(t)),  lat = asin(sin i sin u).
%   toe = toc = refSow rounded to the 16-s resolution of the message.

if nargin < 3 || isempty(sats), sats = navicArtificialDefault(); end
GM     = 3.986005e14;          % as getSatPos.m
OmegaE = 7.2921151467e-5;      % [rad/s], as getSatPos.m
sqrtA  = (GM/OmegaE^2)^(1/6);
toe    = 16*round(refSow/16);

% week number and field layout of a real record
real = find(arrayfun(@(e) ~isempty(e.sqrta) && e.sqrta > 0, eph), 1);
if isempty(real), error('navicArtificialEph: no real NavIC ephemeris to take the week number from.'); end
week = eph(real).weekNrm;

for k = 1:numel(sats)
    s = sats(k);
    if s.prn <= numel(eph) && ~isempty(eph(s.prn).sqrta) && eph(s.prn).sqrta > 0
        error('navicArtificialEph: PRN %d already has a broadcast ephemeris.', s.prn);
    end
    u0 = s.u0Deg*pi/180;
    % Earth-fixed node longitude at toe = lonCross - u0; Omega_k = Omega0 - OmegaE*(tk + toe)
    Omega0 = wrapPi(s.lonCrossDeg*pi/180 - u0 + OmegaE*toe);
    e = eph(real);                                   % copy the field layout
    e.PRN = s.prn;  e.toc = toe;  e.toe = toe;
    e.af0 = 0;  e.af1 = 0;  e.af2 = 0;  e.tgd = 0;
    e.deltan = 0;  e.M0 = wrapPi(u0);  e.ecc = 0;  e.sqrta = sqrtA;
    e.Omega0 = Omega0;  e.omega = 0;  e.Omegadot = 0;
    e.i0 = s.inclDeg*pi/180;  e.idot = 0;
    e.cuc = 0;  e.cus = 0;  e.cic = 0;  e.cis = 0;  e.crc = 0;  e.crs = 0;
    e.svaccur = 0;  e.svhealth = 0;  e.IODEC = 200 + s.prn;  e.weekNrm = week;
    eph(s.prn) = e;
    fprintf('Artificial NavIC PRN %d (%s): i = %.2f deg, crossing %.2f deg E, u0 = %.0f deg\n', ...
        s.prn, s.name, s.inclDeg, s.lonCrossDeg, s.u0Deg);
end
end

function a = wrapPi(a)
a = mod(a + pi, 2*pi) - pi;
end
