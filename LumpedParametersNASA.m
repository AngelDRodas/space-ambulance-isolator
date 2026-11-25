%% LumpedParametersNASA.m  (from LumpedParametersExperimental.m)
% Lumped mass, CG location, and pitch inertia of the
% Patient + Litter (stretcher) + Rack/SMEED + Equipment system,
% followed by a to-scale diagram with labeled lengths, CGs, and inertias.
%
% COORDINATES (from the sketch)
%   z : along the litter, origin at the head end of the litter
%   x : vertical (up), origin at the bottom of the litter
%   y : lateral (into page). J_yy is about the y-axis through the system CG.
%
% Revision 2026-11-19: patient data from NASA-STD-3001 Vol 2 App E (min and
% max crewmember, no regressions); masses from medevac_medical_kits_rev5.xlsx;
% litter length 205 cm; runs both crew cases in one go.
%
% PATIENT DATA (NASA-STD-3001 Vol 2 Appendix E)
%   mass        Table E.5-1        42.64 / 110.22 kg
%   stature     Table E.2-1        148.6 / 194.6 cm (+3% spinal elongation option)
%   J_pat       Table E.5-5, Yp     6.12 / 16.43 kg*m^2 (pitch of a supine patient)
%   CG          Table E.5-3, L(Za) -3.81 / +8.15 cm from the whole-body
%               anatomical origin (+Z toward the head)
%   Whole-body anatomical origin = pelvic axis origin (McConville 1980,
%   AFAMRL-TR-80-119, p.15): midway between the anterior superior iliac
%   spines (ASIS). McConville mean ASIS height = 997.9/1774.9 = 0.562 S,
%   so vertex to ASIS = 0.438 S. Hence
%     z_rel = 0.438 S - L(Za)                    CG below top of head
%   x_rel = CG height above the back. Table E.5-3 L(Xa) = -15.27 / -6.40 cm
%   is measured from the ASIS frontal plane (+X toward the front), so
%     x_rel = buttock depth + L(Xa)
%   Buttock depth scaled from McConville 1980 means: 240.8/1774.9 = 0.1357 S
%   -> 20.2 / 26.4 cm, so x_rel = 4.9 / 20.0 cm. (Check: McConville mean male
%   24.1 - 8.69 = 15.4 cm.) Caveats: NASA's Xa range pairs the female low and
%   male high extremes; soft tissue over the ASIS and buttock compression
%   when supine are ignored.
%   Drawing only: torso thickness = waist depth (15.0 / 30.0 cm), thigh =
%   thigh clearance (13.0 / 20.1 cm, Table E.2-1), calf = 0.69 x thigh
%   (McConville calf/midthigh depth), knee at 0.715 S below the vertex.
%   z_torso = 0.541 S   vertex to gluteal furrow (McConville 1980 mean ratio,
%                       815.1/1774.9 = 0.459 S above floor); the bridge starts here.
%   With spinal elongation on, the 3% S is added above the pelvis, so
%   vertex-to-ASIS and vertex-to-furrow distances grow by 0.03 S.
%
% EQUIPMENT (rev5 Budgets, isolated-mass block, 60.25 kg without patient)
%   stretcher (EMS-SOL vacuum mattress 8.2 + Junkin JSA-200 basket 14.1)  22.3 kg
%   bridge (modified SMEED, assumed)                                      10.0 kg
%   rack devices on the bridge                                            19.5 kg
%   on-patient items (T 1.3 + SC/SP 5.55 + disposables 1.6), added to the
%   patient lump as a point mass at the patient CG                         8.45 kg
%
% COMPONENT CGs
%   z_stretcher = L_s/2                     x_stretcher = H_s/2
%   z_patient   = z_HO + z_rel              x_patient   = H_s + x_rel
%   z_rack      = z_HO + z_torso + L_r/2    x_rack      = H_s + H_r/2
%   z_equip     = z_HO + z_torso + L_e/2    x_equip     = H_s + H_r + H_e/2
%
%   I_rect = (1/12) m (b^2 + h^2)
%   J_yy   = sum( I_i + m_i [ (x_i - x_cg)^2 + (z_i - z_cg)^2 ] )
%
% Requires MATLAB R2016b or newer (local functions in a script).

clear; clc; close all;

%% ---- Inputs ------------------------------------------------------------
% Patient cases (NASA-STD-3001 Vol 2 App E): min = small female, max = large male
P(1).name = 'Min crewmember';  P(1).m_kg = 42.64;   P(1).S_cm = 148.6;
P(1).J_kgm2 = 6.12;   P(1).Za_cm = -3.81;  P(1).Xa_cm = -15.27;
P(1).waistDepth_cm = 15.0;  P(1).thigh_cm = 13.0;
P(2).name = 'Max crewmember';  P(2).m_kg = 110.22;  P(2).S_cm = 194.6;
P(2).J_kgm2 = 16.43;  P(2).Za_cm =  8.15;  P(2).Xa_cm =  -6.40;
P(2).waistDepth_cm = 30.0;  P(2).thigh_cm = 20.1;

p.spinalElongation = true;   % add 3% of stature (microgravity), above the pelvis
p.ASIS_ratio   = 0.562;      % ASIS height / stature (McConville 1980 means)
p.furrow_ratio = 0.541;      % vertex-to-gluteal-furrow / stature (McConville 1980)
p.buttock_ratio = 0.1357;    % buttock depth / stature (McConville 1980 means)
p.knee_ratio   = 0.715;      % vertex-to-knee / stature (drawing only)

% Design parameters [in]
p.stretcherLength_in = 205/2.54;  % 205 cm (Angel 2026-11-19; final size: further study)
p.rackHeight_in      = 8.5;  % NASA hip clearance with margin
p.rackLength_in      = 18;
p.equipHeight_in     = 10.4;     % Propaq height
p.equipLength_in     = 18;       % set by SMEED
p.headrestOffset_in  = 0.5;
p.stretcherHeight_in = 0.75;     % litter thickness in x (kept from previous version)

% Masses [kg] (rev5 Budgets)
p.stretcherMass_kg = 22.3;   % vacuum mattress + basket
p.rackMass_kg      = 10.0;   % bridge (assumed)
p.equipMass_kg     = 19.5;   % rack devices on the bridge
p.onPatientMass_kg = 8.45;   % added to the patient lump

makePlot = true;
checkParams(p);

%% ---- Unit conversions --------------------------------------------------
in2m = 0.0254;  lb2kg = 0.45359237;  cm2m = 0.01;

%% ---- Design geometry [m] -----------------------------------------------
Ls  = p.stretcherLength_in * in2m;   % along z
Hs  = p.stretcherHeight_in * in2m;   % along x (vertical thickness)
Lr  = p.rackLength_in      * in2m;
Hr  = p.rackHeight_in      * in2m;
Le  = p.equipLength_in     * in2m;
He  = p.equipHeight_in     * in2m;
zHO = p.headrestOffset_in  * in2m;

for ic = 1:numel(P)
%% ---- Patient -----------------------------------------------------------
S_cm  = P(ic).S_cm;
dS_cm = 0.03*S_cm*p.spinalElongation;          % elongation above the pelvis
S_eff = S_cm + dS_cm;                           % stature on the litter [cm]
if zHO + S_eff*cm2m > Ls
    warning('%s (%.1f cm) does not fit on a %.1f cm litter.', P(ic).name, S_eff, Ls/cm2m);
end
m_pat    = P(ic).m_kg + p.onPatientMass_kg;     % patient lump incl. on-patient items
zRel_pat = ((1 - p.ASIS_ratio)*S_cm + dS_cm - P(ic).Za_cm) * cm2m;
xRel_pat = (p.buttock_ratio*S_cm + P(ic).Xa_cm) * cm2m;   % CG above the back
zTorso   = (p.furrow_ratio*S_cm + dS_cm) * cm2m;
I_pat    = P(ic).J_kgm2;   % NASA Yp; on-patient items treated as a point mass at the CG

%% ---- Other component masses [kg] ---------------------------------------
m_str = p.stretcherMass_kg;
m_rck = p.rackMass_kg;
m_eqp = p.equipMass_kg;

%% ---- Component table ---------------------------------------------------
name = {'Patient'; 'Stretcher'; 'Bridge/SMEED'; 'Equipment'};

m = [m_pat; m_str; m_rck; m_eqp];

z = [ zHO + zRel_pat;          % patient
      Ls/2;                    % stretcher
      zHO + zTorso + Lr/2;     % rack/SMEED  (starts at the hips)
      zHO + zTorso + Le/2 ];   % equipment   (starts at the hips)

x = [ Hs + xRel_pat;           % patient lies on top of stretcher
      Hs/2;                    % stretcher
      Hs + Hr/2;               % rack/SMEED
      Hs + Hr + He/2 ];        % equipment on top of rack

I = [ I_pat;
      rectI(m_str, Ls, Hs);
      rectI(m_rck, Lr, Hr);
      rectI(m_eqp, Le, He) ];

%% ---- System properties -------------------------------------------------
M    = sum(m);
z_cg = sum(m .* z) / M;
x_cg = sum(m .* x) / M;

d2   = (x - x_cg).^2 + (z - z_cg).^2;
md2  = m .* d2;
J    = sum(I + md2);

out(ic).case = P(ic).name;
out(ic).m_total = M;
out(ic).z_cg    = z_cg;
out(ic).x_cg    = x_cg;
out(ic).J_yy    = J;
out(ic).components = struct('name', name, 'm_kg', num2cell(m), ...
    'z_m', num2cell(z), 'x_m', num2cell(x), ...
    'Icg_kgm2', num2cell(I), 'md2_kgm2', num2cell(md2));

%% ---- Print summary -----------------------------------------------------
fprintf('\n%s: S = %.1f cm (%.1f cm on litter), body mass %.2f kg + %.2f kg on-patient items\n', ...
        P(ic).name, S_cm, S_eff, P(ic).m_kg, p.onPatientMass_kg);
fprintf('%-11s %9s %9s %9s %12s %12s\n', 'Component', 'm [kg]', ...
        'z [m]', 'x [m]', 'I_cg[kgm2]', 'm*d^2[kgm2]');
for k = 1:numel(m)
    fprintf('%-11s %9.3f %9.4f %9.4f %12.4f %12.4f\n', ...
            name{k}, m(k), z(k), x(k), I(k), md2(k));
end
fprintf('\nSYSTEM\n');
fprintf('  Mass        = %8.3f kg   (%8.2f lb)\n', M, M/lb2kg);
fprintf('  z_cg        = %8.4f m    (%8.2f in from head end)\n', z_cg, z_cg/in2m);
fprintf('  x_cg        = %8.4f m    (%8.2f in above litter bottom)\n', x_cg, x_cg/in2m);
fprintf('  J_yy (@CG)  = %8.4f kg*m^2 (%8.1f lb*in^2)\n\n', J, J/(lb2kg*in2m^2));

%% ========================================================================
%  Diagram (all SI: lengths in m, mass in kg, inertia in kg*m^2)
%  ========================================================================
if makePlot
    zc = z';   xc = x';   Ic = I';   mc = m';   md2c = md2';   % already SI
    zS = z_cg; xS = x_cg;

    S_m     = S_eff * cm2m;
    headLen = 0.0938 * S_m;             % head length (McConville ratio)
    torsoThk = P(ic).waistDepth_cm * cm2m;      % drawing only
    thighThk = P(ic).thigh_cm      * cm2m;
    calfThk  = 0.69 * thighThk;
    zKnee    = zHO + (p.knee_ratio*S_cm + dS_cm) * cm2m;
    if Hs + thighThk > Hs + Hr
        warning('%s: thigh (%.3f m) is taller than the bridge clearance (%.3f m).', ...
                P(ic).name, thighThk, Hr);
    end
    topX    = Hs + Hr + He;             % top of equipment

    figure('Color', 'w', 'Position', [80 80 1300 700], 'Name', P(ic).name);
    ax = axes; hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');
    set(ax, 'FontSize', 9, 'Layer', 'top');

    colDim = [0.25 0.25 0.25];
    colSys = [0.85 0.10 0.10];
    colLdr = [0.60 0.60 0.60];

    % Stretcher
    rectangle('Position', [0 0 Ls Hs], 'FaceColor', [0.85 0.85 0.85], ...
              'EdgeColor', 'k', 'LineWidth', 1.5);

    % Patient (head + body)
    rectangle('Position', [zHO, Hs, headLen, headLen], 'Curvature', [1 1], ...
              'FaceColor', [0.80 0.88 1.00], 'EdgeColor', [0.2 0.3 0.6], 'LineWidth', 1.2);
    pcol = {'FaceColor', [0.80 0.88 1.00], 'EdgeColor', [0.2 0.3 0.6], 'LineWidth', 1.2};
    zFur = zHO + zTorso;                                   % gluteal furrow
    rectangle('Position', [zHO + headLen, Hs, zFur - zHO - headLen, torsoThk], ...
              'Curvature', [0.2 0.5], pcol{:});            % torso
    rectangle('Position', [zFur - 0.02, Hs, zKnee - zFur + 0.04, thighThk], ...
              'Curvature', [0.3 0.6], pcol{:});            % thighs (pass under the bridge)
    rectangle('Position', [zKnee, Hs, zHO + S_m - zKnee, calfThk], ...
              'Curvature', [0.3 0.6], pcol{:});            % lower legs

    % Rack / SMEED (open bridge frame, so the patient shows through)
    rectangle('Position', [zc(3) - Lr/2, Hs, Lr, Hr], ...
              'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 2);

    % Equipment
    rectangle('Position', [zc(4) - Le/2, Hs + Hr, Le, He], ...
              'FaceColor', [1.00 0.90 0.70], 'EdgeColor', 'k', 'LineWidth', 1.5);

    % ---- Dimensions [m] ----
    dimH(0, Ls,    -0.10, sprintf('Stretcher length = %.3f m', Ls), colDim);
    dimH(0, zc(1), -0.20, sprintf('z_{patient,CG} = z_{HO} + z_{rel} = %.3f m', zc(1)), colDim);
    dimH(0, zc(3), -0.30, sprintf('z_{SMEED,CG} = z_{HO} + z_{torso} + L_{rack}/2 = %.3f + %.3f + %.3f = %.3f m', ...
                                  zHO, zTorso, Lr/2, zc(3)), colDim);
    dimH(0, zS,    -0.40, sprintf('z_{sys,CG} = %.3f m', zS), colSys);
    dimH(zHO, zHO + zTorso, -0.50, ...
         sprintf('z_{torso} = vertex to gluteal furrow = %.3f m', zTorso), colDim);

    dimV(-0.075, 0, Hs, sprintf('%.3f m', Hs), colDim, 'left');
    dimV(-0.23,  0, xS, sprintf('x_{sys,CG} = %.3f m', xS), colSys, 'left');

    zDimR = zc(3) + Lr/2 + 0.05;
    dimV(zDimR, Hs,      Hs + Hr, sprintf('%.3f m', Hr), colDim, 'right');
    dimV(zDimR, Hs + Hr, topX,    sprintf('%.3f m', He), colDim, 'right');

    dimH(zc(4) - Le/2, zc(4) + Le/2, topX + 0.05, sprintf('%.3f m', Le), colDim);

    text(zHO, Hs + headLen + 0.025, sprintf('z_{HO} = %.3f m', zHO), ...
         'FontSize', 8, 'Color', colDim, 'VerticalAlignment', 'bottom');

    % ---- Callouts (mass, CG, inertia) ----
    zCol = Ls + 0.15;                      % right-hand callout column
    lab  = [ 0.10,  topX + 0.13;           % patient   (upper left)
             zCol,  0.025;                 % stretcher
             zCol,  0.34;                  % rack
             zCol,  0.66 ];                % equipment

    for k = 1:4
        plot([zc(k) lab(k,1)], [xc(k) lab(k,2)], '-', 'Color', colLdr);
        txt = { ['\bf' name{k}], ...
                sprintf('m = %.1f kg', mc(k)), ...
                sprintf('CG (z, x) = (%.3f, %.3f) m', zc(k), xc(k)), ...
                sprintf('I_{cg} = %.3f kg\\cdotm^2', Ic(k)), ...
                sprintf('m d^2 = %.3f kg\\cdotm^2', md2c(k)) };
        text(lab(k,1), lab(k,2), txt, 'FontSize', 8.5, 'BackgroundColor', 'w', ...
             'EdgeColor', colLdr, 'Margin', 2, 'VerticalAlignment', 'middle');
    end

    for k = 1:4
        cgMarker(zc(k), xc(k), 0.025, 'k');
    end
    cgMarker(zS, xS, 0.033, colSys);

    % ---- System summary ----
    sumTxt = { '\bfSYSTEM (about y-axis through system CG)', ...
               sprintf('Patient: S = %.3f m, m = %.1f kg (incl. %.2f kg items)', S_m, m_pat, p.onPatientMass_kg), ...
               sprintf('M_{total} = %.1f kg', M), ...
               sprintf('CG (z, x) = (%.3f, %.3f) m', zS, xS), ...
               sprintf('J_{yy} = \\Sigma(I_{cg} + m d^2) = %.3f kg\\cdotm^2', J) };
    text(zCol, -0.38, sumTxt, 'FontSize', 9, 'Color', colSys, ...
         'BackgroundColor', 'w', 'EdgeColor', colSys, 'Margin', 3, ...
         'VerticalAlignment', 'middle');

    % ---- Coordinate triad ----
    z0 = -0.50; x0 = topX + 0.10; aL = 0.13;
    plot([z0 z0+aL], [x0 x0], 'k-', 'LineWidth', 1.2);
    plot(z0+aL, x0, 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    plot([z0 z0], [x0 x0+aL], 'k-', 'LineWidth', 1.2);
    plot(z0, x0+aL, 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    text(z0+aL+0.025, x0, 'z', 'FontSize', 10);
    text(z0, x0+aL+0.04, 'x', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(z0+0.04, x0+0.05, '\otimes y', 'FontSize', 9);

    xlim([-0.61, Ls + 0.86]);
    ylim([-0.58, topX + 0.36]);
    xlabel('z [m]  (from head end of litter)');
    ylabel('x [m]  (from bottom of litter)');
    title(sprintf('Lumped mass properties (%s): Patient + Litter + Bridge + Equipment', P(ic).name));
    grid on; set(ax, 'GridAlpha', 0.12);
end
end   % patient cases

%% ---- Summary over cases (for the pitch-heave model) ---------------------
fprintf('\n%-16s %9s %9s %9s %12s\n', 'Case', 'M [kg]', 'z_cg [m]', 'x_cg [m]', 'J_yy [kgm2]');
for ic = 1:numel(out)
    fprintf('%-16s %9.2f %9.4f %9.4f %12.3f\n', out(ic).case, out(ic).m_total, ...
            out(ic).z_cg, out(ic).x_cg, out(ic).J_yy);
end

%% ========================================================================
%  Local functions
%  ========================================================================
function I = rectI(m, base, height)
% Centroidal moment of inertia of a uniform rectangle: (1/12) m (b^2 + h^2)
I = m * (base^2 + height^2) / 12;
end

function checkParams(p)
need = {'stretcherMass_kg', 'rackMass_kg', 'equipMass_kg', 'onPatientMass_kg', 'stretcherHeight_in'};
for k = 1:numel(need)
    v = p.(need{k});
    if ~isscalar(v) || isnan(v)
        error('Parameter p.%s is not defined yet.', need{k});
    end
end
end

function dimH(z1, z2, x, label, col)
% Horizontal dimension line with arrowheads and centered label (units: m)
plot([z1 z2], [x x], '-', 'Color', col, 'LineWidth', 1);
plot(z1, x, '<', 'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 4);
plot(z2, x, '>', 'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 4);
plot([z1 z1], [x-0.02 x+0.02], '-', 'Color', col);
plot([z2 z2], [x-0.02 x+0.02], '-', 'Color', col);
text((z1+z2)/2, x + 0.01, label, 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', col);
end

function dimV(z, x1, x2, label, col, side)
% Vertical dimension line with arrowheads; label on the given side (units: m)
plot([z z], [x1 x2], '-', 'Color', col, 'LineWidth', 1);
plot(z, x1, 'v', 'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 4);
plot(z, x2, '^', 'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 4);
if strcmp(side, 'left')
    text(z - 0.015, (x1+x2)/2, label, 'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', col);
else
    text(z + 0.015, (x1+x2)/2, label, 'HorizontalAlignment', 'left', 'FontSize', 8, 'Color', col);
end
end

function cgMarker(zc, xc, r, col)
% Standard CG symbol: circle with two filled opposite quadrants
t = linspace(0, 2*pi, 60);
patch(zc + r*cos(t), xc + r*sin(t), 'w', 'EdgeColor', col, 'LineWidth', 1.2);
t1 = linspace(0,      pi/2, 20);
t3 = linspace(pi, 3*pi/2,   20);
patch([zc, zc + r*cos(t1)], [xc, xc + r*sin(t1)], col, 'EdgeColor', col);
patch([zc, zc + r*cos(t3)], [xc, xc + r*sin(t3)], col, 'EdgeColor', col);
end