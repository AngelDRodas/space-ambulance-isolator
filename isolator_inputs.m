function P = isolator_inputs()
%ISOLATOR_INPUTS  Shared inputs for the stretcher landing isolator study.
% Both model scripts read this file, so the vibrations (filter) model and
% the bond-graph (RK4) model always see the same case. Change inputs here,
% not inside the model scripts.
%
%   ShockAttentuationSpaceAmbulanceBrinkley.m    vibrations EOM, exact discrete filter()
%   ShockAttentuationSpaceAmbulanceBondGraphs.m  bond-graph state equations, RK4
%   compare_models.m                             error comparison of the two
%
% Revision 2026-11-19: isolated mass, patient mass and touchdown velocity
% updated (see notes beside each input).

P.g = 9.81;                      % [m/s^2]


%% ---------------- Isolated mass ----------------
% Everything rigidly fastened to the stretcher that rides on the isolator
% (medevac_medical_kits_rev5.xlsx, Budgets, isolated-mass block):
%   medical equipment on the stretcher   27.95 kg  (rack devices 19.5, on-patient
%                                                   T items 1.3, SC/SP items 5.55,
%                                                   in-use disposables 1.6)
%   vacuum mattress EMS-SOL               8.2  kg
%   basket stretcher Junkin JSA-200      14.1  kg
%   equipment bridge (SMEED-based)       10.0  kg  (assumed)
%   -----------------------------------------------
%   package without patient              60.25 kg
% NOT isolated: E cylinder (capsule mount in flight), LUCAS (station only),
% flotation collar (stowed), main O2 cylinders (vehicle structure mounts).
% Patient: NASA-STD-3001 Vol 2 App E, Table E.5-1, unsuited body mass
% 42.64 kg (lightest) to 110.22 kg (heaviest).
P.m_package = 60.25;                       % [kg]
P.m_patient = [42.64 110.22];              % [kg]
P.m_list    = P.m_package + P.m_patient;   % [kg] -> 102.89 and 170.47
% The light case drives acceleration and Brinkley (a stiff isolator gives the
% light patient a harder ride); the heavy case drives stroke.

%% ---------------- Impact input ----------------
P.Vi      = 8.9;          % touchdown velocity [m/s]: Orion splashdown <= 20 mph
                          % (NASA 2022, Orion reference guide). Was 7.62 m/s (Apollo 25 ft/s).
P.pulse   = 'halfsine';   % 'halfsine' | 'rect' | 'data'
% Peak capsule accelerometer reading in the pulse [g] (synthetic only).
% Reference: Apollo splashdown, ~15 g peak at 7.62 m/s (25 ft/s).
% Scaled to the Orion touchdown velocity assuming the hydrodynamic water-entry
% load goes as rho*V^2 (same heat-shield shape class, similar mass per area):
%   a_pk = 1 + (a_ref - 1)*(Vi/V_ref)^2   -> 20.1 g at 8.9 m/s, half-sine 75 ms.
% Linear scaling (same pulse duration) would give 17.4 g; V^2 is the more
% conservative and more physical choice. Replace with a measured pulse
% (pulse = 'data') when one is available.
P.a_ref   = 15;           % [g]   Apollo reference peak
P.V_ref   = 7.62;         % [m/s] Apollo reference touchdown velocity
P.a_in_pk = 1 + (P.a_ref - 1)*(P.Vi/P.V_ref)^2;   % [g]

% Measured accelerometer data (used when pulse = 'data')
P.data_file  = 'accel.csv';  % two columns: time [s], vertical acceleration
P.data_units = 'g';          % 'g' | 'm/s^2'
P.data_sign  = +1;           % -1 if the pre-impact reading is about -1 g
P.t_imp_data = [];           % impact time in the file's time base; [] = auto-detect
P.detect_g   = 2;            % auto-detect threshold [g]
P.pad_post   = true;         % pad the end with the floating 1 g reading

%% ---------------- Constraints ----------------
P.stroke_max  = 0.50;     % [m]
P.stroke_def  = 'total';  % 'total' | 'dynamic'
% Legacy limits (10 g peak, 2 g*s impulse). The thesis now treats them as
% superseded by Brinkley. Set use_legacy = false to size on stroke + beta only.
P.use_legacy  = false;    % Angel 2026-11-19: size on stroke + Brinkley only
P.accel_max   = 10;       % [g]
P.impulse_max = 2;        % [g*s]
P.accel_def   = 'reading';   % 'reading' | 'kinematic'

%% ---------------- Brinkley Dynamic Response ----------------
P.use_brinkley   = true;
P.brinkley_axis  = 'x';     % supine patient: splashdown load is chest-to-back (X)
P.crew_condition = 'ill';   % 'healthy' | 'ill' (ill/injured limits, ~75% of healthy)
P.brinkley_sign  = +1;      % +1: upward stretcher accel pushes the back -> eyeballs-in (+X)
P.brinkley_1g    = true;    % true: DR = 1 at rest (1 g kept in the input, conservative)
% ACRV landing impact table (Brinkley et al.)
%            wn [rad/s]  zeta   healthy [+ -]   ill/injured [+ -]
P.BTAB.x = {62.8, 0.200, [35.0 -28.0], [26.3 -21.0]};
P.BTAB.y = {58.0, 0.090, [14.0 -14.0], [10.5 -10.5]};
P.BTAB.z = {52.9, 0.224, [15.2  -9.0], [11.4  -6.8]};

%% ---------------- Margins and selection ----------------
P.stroke_margin  = 0.00;
P.accel_margin   = 0.00;
P.impulse_margin = 0.00;
P.beta_margin    = 0.00;
P.selection      = 'stroke_first';   % 'stroke_first' | 'balanced'

%% ---------------- Sweep ----------------
P.k_vec = linspace(1e2, 4e4, 200);   % [N/m]
P.c_vec = linspace(1e2, 4e3, 200);   % [N*s/m]
P.axis_scale  = 'linear';
P.zeta_guides = [0.1 0.3 0.7 1.0];

% Check design: both scripts save its full time history for compare_models.m
% (Design B, sized in October for 100/150 kg at 7.62 m/s).
P.k_check = 13000;   P.c_check = 1700;

%% ---------------- Time grid ----------------
P.t_pre  = 0.2;     % steady descent before impact [s]
P.t_post = 2.5;     % after the pulse [s]
P.n_min  = 200;     % minimum samples across the pulse
P.ppc    = 50;      % samples per period of the fastest system (isolator OR Brinkley)
% Convergence check: double n_min and ppc, rerun both scripts and compare_models.
% If the model-to-model error drops ~4x, the difference is numerical (2nd order).

%% ---------------- Pitch-heave model (PitchHeaveIsolator.m) ----------------
% From LumpedParametersNASA.m (2026-11-19): 205 cm litter, NASA-STD-3001 min/max
% crew, rev5 masses; order matches m_list (min, max).
P.ph_J      = [18.813 32.641];   % [kg m^2] pitch inertia about the system CG
P.ph_zcg    = [0.9052 0.9542];   % [m] system CG from the head end
P.ph_zpat   = [0.7462 0.8419];   % [m] patient CG from the head end (Brinkley point)
P.ph_z_head = 0;                 % [m] head isolator
P.ph_z_foot = 2.05;              % [m] foot isolator (litter end)
P.ph_pulses = {'halfsine'};   % Angel 2026-11-19: half-sine only. Add 'rect' for a pulse-shape sensitivity run (worst case over shapes).
P.ph_t_post = 1.5;               % [s] shorter window keeps the 4-D search fast
P.ph_k_equal = linspace(1000, 30000, 100);   % [N/m]   per isolator
P.ph_c_equal = linspace(50, 3000, 100);      % [N*s/m] per isolator
P.ph_sep_span = [0.5 0.5];       % separate search: +/- fraction around the equal optimum (k, c)
P.ph_sep_n    = [9 9];           % points per axis -> 9^4 = 6561 designs
P.ph_polish   = true;            % fminsearch polish of the separate design
P.ph_ode45_check = true;         % independent Newton-Euler/ode45 check of the pitch coupling
end
