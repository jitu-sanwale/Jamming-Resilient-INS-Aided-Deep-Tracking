%INIT  Jamming-resilient INS-aided deep tracking of NavIC signals (IEEE TAES R1).
%
%   Opens the mode menu (aj\ajMenu.m): jammer-free reference, navigation-level
%   campaign, the two ablations, tables and statistics, jammer localization,
%   array-level figures, or one seed end to end. Each mode prompts for its
%   options. For non-interactive use define a struct AJ before running init,
%   e.g.  AJ = struct('mode', 2, 'worker', 1, 'nWorkers', 2, 'seeds', 1:10);  init
%   (fields: see aj\ajMenu.m).
%
%   Simulator: eight-element (7+1) CRPA IF generation for NavIC L5 SPS on an
%   F-16 (nonlinear 6-DOF model), with a navigation-grade IMU and jammers.
%   Array processing: PI and INS-steered per-satellite MVDR with the body-rate
%   block schedule and blanking. Receiver: INS-aided ultratight NavIC tracking
%   with array aiding. Requirements: MATLAB R2023a or later with the Signal
%   Processing, Communications and Curve Fitting Toolboxes; Python 3 (NumPy,
%   SciPy, Matplotlib) for the localization and the array-level figures.
%
%   The receiver builds on GPSSDR_vt (B. Xu and L.-T. Hsu, GPS Solutions 23,
%   46, 2019) and the simulator on GPS_L1CA_signal_generator, the open-source
%   IF signal simulator of Y. Li (Beijing Information Science and Technology
%   University) as modified by J. Sun; both include routines from SoftGNSS
%   (K. Borre, D. Akos, D. Plausinaitis et al.).
%   The F-16 aerodynamic data are from F16Model.jl (MIT licence, see
%   Simulator\F16Aircraft). The broadcast ephemerides
%   (Simulator\IITK00IND_R_20250011000_01H_MN.rnx) are data of the
%   International GNSS Service (IGS), https://igs.org.
%   The original notices are kept in the file headers. Licence: GPL-3.0-or-later.

pkg = fileparts(mfilename('fullpath'));
addpath(fullfile(pkg, 'aj'));
if exist('AJ', 'var') && isstruct(AJ), ajMenu(AJ); else, ajMenu; end
