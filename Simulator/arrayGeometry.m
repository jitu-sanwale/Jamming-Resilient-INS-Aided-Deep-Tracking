function A = arrayGeometry(opt)
%ARRAYGEOMETRY 7+1 CRPA for NavIC L5: element positions, coupling, steering.
%   A = arrayGeometry(opt)
%
%   Shared by the multi-element IF generator and the array processor, so both
%   use one convention. Elements: 0 at the array centre (reference), 1..7 on a
%   ring of radius lambda_L5/2 in the body x-y plane, element 1 on the body
%   x axis (nose), counted towards body +y (right wing). Body axes as
%   applyAttitudeToBodyFrame.m: x forward, y right, z down; the identical,
%   zenith-pointing patches look along body -z. The array centre is the
%   trajectory point (no lever arm).
%
%   opt (all optional):
%     fc            carrier [Hz]                (default 1176.45e6, NavIC L5)
%     radius        ring radius [m]             (default lambda/2)
%     couplingDb    adjacent-element coupling [dB]; -Inf = none (default -20)
%   Fields of A:
%     M, lambda, fc, posB (3xM, body frame [m]), C (MxM coupling matrix)
%     steer(uB)     M x N steering vectors for N body-frame unit vectors uB (3xN),
%                   including coupling: a = C * exp(-1i*2*pi/lambda*(posB'*uB)).
%                   Sign convention of the IF generator: the replica carrier is
%                   exp(-1i*w*Tx), and an element displaced by d towards the
%                   satellite receives the content of Tx + (u'd)/c.
%
%   Mutual coupling (illustrative, as in the paper): c_mn = kappa (d_adj/d_mn)^2
%   exp(-1i (k d_mn + pi/2)), kappa = 10^(couplingDb/20). The coupling acts on
%   the incident field, i.e. on satellites and jammers, not on receiver noise.

if nargin < 1, opt = struct(); end
c = 299792458;
A.fc = getf(opt, 'fc', 1176.45e6);
A.lambda = c/A.fc;
r = getf(opt, 'radius', A.lambda/2);
A.M = 8;
k = 0:6;
A.posB = [zeros(3,1), [r*cos(2*pi*k/7); r*sin(2*pi*k/7); zeros(1,7)]];
cdb = getf(opt, 'couplingDb', -20);
D = sqrt(sum((permute(A.posB, [2 3 1]) - permute(A.posB, [3 2 1])).^2, 3));
A.C = eye(A.M);
if ~isempty(cdb) && isfinite(cdb)
    kap = 10^(cdb/20);  dAdj = 2*r*sin(pi/7);  K = 2*pi/A.lambda;
    off = ~eye(A.M);
    A.C(off) = kap*(dAdj./D(off)).^2 .* exp(-1i*(K*D(off) + pi/2));
end
A.couplingDb = cdb;
lam = A.lambda;  P = A.posB;  Cm = A.C;
A.steer = @(uB) Cm*exp(-1i*2*pi/lam*(P'*uB));
end

function v = getf(S, f, d)
if isfield(S, f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end
