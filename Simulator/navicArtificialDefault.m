function sats = navicArtificialDefault()
%NAVICARTIFICIALDEFAULT Three artificial satellites that complete the nominal
%   seven-satellite NavIC constellation (3 GEO + 4 IGSO) around the four
%   satellites of the broadcast ephemeris file (PRN 2, 6, 9 IGSO/GEO at the
%   55 deg E and 32.5/129.5 deg E slots, PRN 10): the vacant slots are the GEO
%   at 83 deg E and the IGSO pair crossing the equator at 111.75 deg E
%   (inclination 29 deg), placed on the northern and southern lobes.
%   PRN 12-14 avoid the PRNs of real NavIC satellites.
sats = struct('name', {'A1 GEO 83E', 'A2 IGSO 111.75E N', 'A3 IGSO 111.75E S'}, ...
              'prn', {12, 13, 14}, ...
              'inclDeg', {0, 29, 29}, ...
              'lonCrossDeg', {83, 111.75, 111.75}, ...
              'u0Deg', {0, 90, 270});
end
