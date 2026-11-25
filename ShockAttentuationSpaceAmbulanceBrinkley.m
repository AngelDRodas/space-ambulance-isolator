%% ShockAttentuationSpaceAmbulanceBrinkley.m  (was splashdown_sweep_brinkley.m)
% Revision 2026-11-19: inputs moved to isolator_inputs.m (shared with the
% bond-graph script); pulse generation moved to isolator_make_input.m;
% saves results_vibrations.mat for compare_models.m.
%
% Parameter sweep of a base-excited SDOF isolated-stretcher model during
% capsule splashdown. Finds the spring (k) and damper (c) constants that
% satisfy stroke, acceleration, impulse and Brinkley limits for each mass
% (isolated package + patient),
% and the design space that satisfies ALL masses at once.
%
%   EOM:  m*zdd + c*zd + k*z = -m*a(t),   z = x - y  (x, y positive UP)
%
% a(t) is the capsule accelerometer reading (specific force): +1 g under
% chutes, a pulse at splashdown, then +1 g floating. Gravity is contained
% in a(t), so the model starts in equilibrium at z(0) = -m*a0/k.
%
% Input: a synthetic pulse ('halfsine' or 'rect') or measured accelerometer
% data ('data'). With data, Vi and a_in_pk are not used; the record itself
% sets the pulse, and the pre-impact average sets the equilibrium.
%
% Synthetic pulse duration is set so the pulse removes the touchdown velocity Vi:
%   half-sine:  T = pi*Vi/(2*A)      rectangle:  T = Vi/A
%   with A = (a_in_pk - 1)*g, the pulse height above the 1 g baseline.
%
% Metrics (edit run_case at the bottom to change any definition):
%   stroke  : 'total'   = max compression from spring free length
%                         (includes the 1 g static preload m*g/k)
%             'dynamic' = max travel beyond the static position
%   accel   : 'reading'   = peak |accelerometer reading on the mass| (1 g at rest)
%             'kinematic' = peak |reading - 1 g|
%   impulse : integral of |reading - 1 g| dt from impact to end of window [g*s]
%             (lower bound is Vi/g; rebound oscillations add to it)
%   beta    : Brinkley Dynamic Response criterion (acceptable impact: beta <= 1)
%
% Brinkley DR model: the occupant's body is a second SDOF system driven by
% the acceleration of the seat (the stroking mass in this model):
%     d'' + 2*zeta_B*wn_B*d' + wn_B^2*d = A(t),    DR(t) = wn_B^2*d(t)/g
% Single axis:  beta = max over t of |DR(t)| / DR_limit (limit picked by sign).
% Parameters and limits: ACRV landing impact table (Brinkley et al.).
%
% No toolboxes required: each case is simulated with an exact discrete
% state-space model applied via filter().


clear; clc; close all

%% ======================= INPUTS (shared) =======================
% All inputs live in isolator_inputs.m so this script and the bond-graph
% script always run the same case. Edit that file, not this one.
IN = isolator_inputs();   % (named IN: P is used below for pass masks)
g = IN.g;
Vi = IN.Vi;  pulse = IN.pulse;  a_in_pk = IN.a_in_pk;
stroke_max = IN.stroke_max;  stroke_def = IN.stroke_def;  accel_def = IN.accel_def;
if IN.use_legacy
    accel_max = IN.accel_max;  impulse_max = IN.impulse_max;
else
    accel_max = Inf;          impulse_max = Inf;     % legacy limits off
end
use_brinkley = IN.use_brinkley;  brinkley_axis = IN.brinkley_axis;
crew_condition = IN.crew_condition;  brinkley_sign = IN.brinkley_sign;
brinkley_1g = IN.brinkley_1g;  BTAB = IN.BTAB;
stroke_margin = IN.stroke_margin;  accel_margin = IN.accel_margin;
impulse_margin = IN.impulse_margin;  beta_margin = IN.beta_margin;
selection = IN.selection;
m_list = IN.m_list;  k_vec = IN.k_vec;  c_vec = IN.c_vec;
axis_scale = IN.axis_scale;  zeta_guides = IN.zeta_guides;
fprintf('Isolated mass: %.2f kg package + patient %s kg -> m = %s kg\n', ...
        IN.m_package, mat2str(IN.m_patient), mat2str(m_list));

%% ======================= INPUT PULSE =======================
[t, a, a0, t_imp, dt, T, post] = isolator_make_input(IN);
u    = a - a0;                  % input deviation from steady descent
u    = [(u(1:end-1) + u(2:end))/2; u(end)];   % interval average: 2nd-order accurate

if strcmp(pulse, 'data')
    fprintf('Data: %s, pre-impact reading %.3f g, peak %.2f g, dt = %.3f ms\n', ...
            IN.data_file, a0/g, max(abs(a))/g, 1e3*dt);
    fprintf('Impact at t = %.4f s; velocity change from data: %.3f m/s\n\n', ...
            t_imp, trapz(t(post), a(post) - a0));
else
    fprintf('Pulse: %s, %.1f g peak, T = %.1f ms, dt = %.3f ms\n', ...
            pulse, a_in_pk, 1e3*T, 1e3*dt);
    fprintf('Velocity removed: %.3f m/s (Vi = %.3f)\n\n', trapz(t, a - a0), Vi);
end

%% ======================= BRINKLEY PARAMETERS =======================
% ACRV landing impact dynamic response limits (table in isolator_inputs.m)
bt = BTAB.(brinkley_axis);
BR.on   = use_brinkley;
BR.wn   = bt{1};  BR.zeta = bt{2};
switch crew_condition
    case 'healthy', BR.lim = bt{3};
    case 'ill',     BR.lim = bt{4};
    otherwise, error('crew_condition must be ''healthy'' or ''ill''');
end
BR.sgn = brinkley_sign;  BR.inc1g = brinkley_1g;
fprintf('Brinkley %s-axis, %s crew: wn = %.1f rad/s, zeta = %.3f, DR limits +%.1f / %.1f\n\n', ...
        upper(brinkley_axis), crew_condition, BR.wn, BR.zeta, BR.lim(1), BR.lim(2));

%% ======================= SWEEP =======================
nm = numel(m_list);  nk = numel(k_vec);  nc = numel(c_vec);
[KK, CC] = meshgrid(k_vec, c_vec);        % rows: c, columns: k

STROKE = zeros(nc, nk, nm);
ACC    = zeros(nc, nk, nm);
IMP    = zeros(nc, nk, nm);
BETA   = zeros(nc, nk, nm);

tic
for im = 1:nm
    m = m_list(im);
    for ic = 1:nc
        for ik = 1:nk
            [STROKE(ic,ik,im), ACC(ic,ik,im), IMP(ic,ik,im), BETA(ic,ik,im)] = ...
                run_case(m, k_vec(ik), c_vec(ic), u, a0, t, post, dt, g, ...
                         stroke_def, accel_def, BR);
        end
    end
end
fprintf('Swept %d designs x %d masses in %.1f s\n\n', nk*nc, nm, toc);

% Brinkley is ignored (beta set to 0) when use_brinkley = false
if ~use_brinkley, BETA(:) = 0; end

% Hard limits
PASS     = STROKE <= stroke_max & ACC <= accel_max & IMP <= impulse_max & BETA <= 1;
PASS_ALL = all(PASS, 3);

% Margin-adjusted limits
s_lim = (1 - stroke_margin)  * stroke_max;
a_lim = (1 - accel_margin)   * accel_max;
i_lim = (1 - impulse_margin) * impulse_max;
b_lim = (1 - beta_margin);
PASS_M     = STROKE <= s_lim & ACC <= a_lim & IMP <= i_lim & BETA <= b_lim;
PASS_M_ALL = all(PASS_M, 3);

% Worst ratio to the margin-adjusted limits: <= 1 means all margins met.
% Used to rank designs that meet every margin, and by the 'balanced' rule.
RATIO = max(cat(4, STROKE/s_lim, ACC/a_lim, IMP/i_lim, BETA/b_lim), [], 4);
WORST = max(RATIO, [], 3);                % worst over all masses

%% ======================= SUMMARY =======================
for im = 1:nm
    P = PASS(:,:,im);
    fprintf('m = %g kg: %d of %d designs meet limits, %d also meet margins\n', ...
            m_list(im), nnz(P), numel(P), nnz(PASS_M(:,:,im)));
    if any(P(:))
        fprintf('   k range: %.0f - %.0f N/m\n',   min(KK(P)), max(KK(P)));
        fprintf('   c range: %.0f - %.0f N*s/m\n', min(CC(P)), max(CC(P)));
    end
end

S_W   = max(STROKE, [], 3);                          % worst stroke over masses
OTHER = max(max(max(ACC/a_lim, IMP/i_lim), BETA/b_lim), [], 3);  % worst accel/impulse/beta ratio
switch selection
    case 'balanced'
        obj = WORST;
        sel_txt = 'balanced: lowest worst ratio to margin-adjusted limits';
    case 'stroke_first'
        if any(PASS_M_ALL(:))
            obj = WORST;  obj(~PASS_M_ALL) = Inf;
            sel_txt = 'meets every margin for every mass';
        elseif any(S_W(:) <= s_lim)
            obj = OTHER;  obj(S_W > s_lim) = Inf;
            sel_txt = 'stroke margin kept; accel/impulse margins traded';
        elseif any(S_W(:) <= stroke_max)
            obj = OTHER;  obj(S_W > stroke_max) = Inf;
            sel_txt = 'stroke limit kept, stroke margin not achievable';
        else
            obj = S_W;
            sel_txt = 'every design bottoms out; least stroke chosen';
        end
    otherwise
        error('selection must be ''stroke_first'' or ''balanced''');
end
[~, ib]  = min(obj(:));
k_best   = KK(ib);  c_best = CC(ib);
fprintf('\nCommon to all masses: %d meet limits, %d also meet margins\n', ...
        nnz(PASS_ALL), nnz(PASS_M_ALL));
fprintf('Selected design (%s):\n  k = %.0f N/m, c = %.0f N*s/m\n', ...
        sel_txt, k_best, c_best);
fprintf('Margin = remaining fraction of each hard limit (negative = exceeded)\n');
fprintf('%8s %7s %10s %8s %9s %8s %12s %8s %7s %8s %8s\n', 'm [kg]', 'zeta', ...
        'stroke[m]', 'margin', 'accel[g]', 'margin', 'impulse[g*s]', 'margin', ...
        'beta', 'margin', 'DR peak');
for im = 1:nm
    [s, ac, ip, bt_, ~, ~, DR] = run_case(m_list(im), k_best, c_best, u, a0, t, post, dt, g, ...
                                          stroke_def, accel_def, BR);
    fprintf('%8g %7.3f %10.3f %7.0f%% %9.2f %7.0f%% %12.3f %7.0f%% %7.3f %7.0f%% %8.2f\n', m_list(im), ...
            c_best/(2*sqrt(k_best*m_list(im))), ...
            s,  100*(1 - s/stroke_max), ac, 100*(1 - ac/accel_max), ...
            ip, 100*(1 - ip/impulse_max), bt_, 100*(1 - bt_), max(DR));
end

%% ======================= CSV OUTPUT =======================
MM3 = repmat(reshape(m_list, 1, 1, nm), nc, nk);
KK3 = repmat(KK, 1, 1, nm);
CC3 = repmat(CC, 1, 1, nm);
ZZ3 = CC3 ./ (2*sqrt(KK3 .* MM3));
results = table(MM3(:), KK3(:), CC3(:), ZZ3(:), STROKE(:), ACC(:), IMP(:), BETA(:), ...
    PASS(:), PASS_M(:), RATIO(:), ...
    'VariableNames', {'mass_kg','k_Npm','c_Nspm','zeta','stroke_m', ...
                      'accel_g','impulse_gs','brinkley_beta','pass_limits','pass_margins','margin_ratio'});
writetable(results, 'sweep_results.csv');

if any(PASS_ALL(:))
    common = table(KK(PASS_ALL), CC(PASS_ALL), PASS_M_ALL(PASS_ALL), WORST(PASS_ALL), ...
        'VariableNames', {'k_Npm','c_Nspm','meets_margins','worst_margin_ratio'});
    common = sortrows(common, 'worst_margin_ratio');
    writetable(common, 'feasible_all_masses.csv');
end

% MAT file for compare_models.m: full sweep plus time histories of the
% check design (same k, c in both models; empty = this model's selection).
if isempty(IN.k_check), IN.k_check = k_best;  IN.c_check = c_best;  end
z_c = zeros(numel(t), nm);  am_c = z_c;  DR_c = z_c;
for im = 1:nm
    [~, ~, ~, ~, z_c(:,im), am_c(:,im), DR_c(:,im)] = run_case(m_list(im), ...
        IN.k_check, IN.c_check, u, a0, t, post, dt, g, stroke_def, accel_def, BR);
end
R = struct('model', 'vibrations (filter)', 'IN', IN, 'm_list', m_list, ...
           'k_vec', k_vec, 'c_vec', c_vec, 'STROKE', STROKE, 'ACC', ACC, ...
           'IMP', IMP, 'BETA', BETA, 'PASS', PASS, 'k_best', k_best, ...
           'c_best', c_best, 't', t, 'a', a, 'dt', dt, ...
           'chk', struct('k', IN.k_check, 'c', IN.c_check, 'stroke', -z_c, ...
                         'am', am_c, 'DR', DR_c));
save('results_vibrations.mat', '-struct', 'R');

%% ======================= PLOTS: EACH MASS =======================
pass_col = [0.72 0.88 0.72];
for im = 1:nm
    m = m_list(im);
    figure('Name', sprintf('m = %g kg', m), 'Color', 'w'); hold on
    P = PASS(:,:,im);

    PM = PASS_M(:,:,im);
    hP = plot([KK(P & ~PM); NaN], [CC(P & ~PM); NaN], 'o', 'MarkerSize', 4, ...
              'MarkerFaceColor', [0.65 0.85 0.65], 'MarkerEdgeColor', 'none');
    hM = plot([KK(PM); NaN], [CC(PM); NaN], 'o', 'MarkerSize', 4, ...
              'MarkerFaceColor', [0.1 0.5 0.1], 'MarkerEdgeColor', 'none');
    if stroke_margin > 0
        contour(KK, CC, STROKE(:,:,im), [s_lim s_lim], '--', ...
                'LineColor', [0.85 0.2 0.2], 'HandleVisibility', 'off');
    end
    if accel_margin > 0
        contour(KK, CC, ACC(:,:,im), [a_lim a_lim], '--', ...
                'LineColor', [0.15 0.35 0.8], 'HandleVisibility', 'off');
    end
    [~, hS] = contour(KK, CC, STROKE(:,:,im), [stroke_max stroke_max], ...
                      'LineColor', [0.85 0.2 0.2], 'LineWidth', 1.6);
    if isfinite(accel_max)
        [~, hA] = contour(KK, CC, ACC(:,:,im), [accel_max accel_max], ...
                          'LineColor', [0.15 0.35 0.8], 'LineWidth', 1.6);
        [~, hI] = contour(KK, CC, IMP(:,:,im), [impulse_max impulse_max], ...
                          'LineColor', [0.9 0.55 0.1], 'LineWidth', 1.6);
    else            % legacy limits off: placeholders keep the legend intact
        hA = plot(NaN, NaN, 'Color', [0.15 0.35 0.8], 'LineWidth', 1.6);
        hI = plot(NaN, NaN, 'Color', [0.9 0.55 0.1], 'LineWidth', 1.6);
    end
    if use_brinkley
        [~, hBr] = contour(KK, CC, BETA(:,:,im), [1 1], ...
                           'LineColor', [0.55 0.2 0.7], 'LineWidth', 1.6);
    else
        hBr = plot(NaN, NaN, 'Color', [0.55 0.2 0.7], 'LineWidth', 1.6);
    end
    add_zeta_lines(k_vec, c_vec, m, zeta_guides);

    [~, ibm] = min(reshape(RATIO(:,:,im), [], 1));
    hB = plot(KK(ibm), CC(ibm), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 13);

    set(gca, 'XScale', axis_scale, 'YScale', axis_scale, 'Layer', 'top'); box on
    xlim([min(k_vec) max(k_vec)]); ylim([min(c_vec) max(c_vec)]);
    xlabel('Spring constant k [N/m]'); ylabel('Damper constant c [N\cdot{}s/m]');
    title(sprintf('m = %g kg: %d meet limits, %d meet margins', m, nnz(P), nnz(PM)));
    legend([hP hM hS hA hI hBr hB], {'Meets limits only', 'Meets limits + margins', ...
        sprintf('Stroke = %.2f m', stroke_max), ...
        sprintf('Accel = %g g', accel_max), ...
        sprintf('Impulse = %g g\\cdot{}s', impulse_max), ...
        sprintf('Brinkley \\beta = 1 (%s, %s)', upper(brinkley_axis), crew_condition), ...
        'Best for this mass'}, 'Location', 'northeast');
end

%% ======================= PLOT: ALL MASSES =======================
figure('Name', 'Common solution space', 'Color', 'w'); hold on
cols = lines(nm);
h = gobjects(0);  lab = {};
LO = PASS_ALL & ~PASS_M_ALL;
h(end+1) = plot([KK(LO); NaN], [CC(LO); NaN], 'o', 'MarkerSize', 4, ...
                'MarkerFaceColor', [0.65 0.85 0.65], 'MarkerEdgeColor', 'none');
lab{end+1} = sprintf('All masses, limits only (%d)', nnz(LO));
h(end+1) = plot([KK(PASS_M_ALL); NaN], [CC(PASS_M_ALL); NaN], 'o', 'MarkerSize', 4, ...
                'MarkerFaceColor', [0.1 0.5 0.1], 'MarkerEdgeColor', 'none');
lab{end+1} = sprintf('All masses, limits + margins (%d)', nnz(PASS_M_ALL));
for im = 1:nm
    if any(any(PASS(:,:,im)))
        [~, hm] = contour(KK, CC, double(PASS(:,:,im)), [0.5 0.5], ...
                          'LineColor', cols(im,:), 'LineWidth', 1.8);
        h(end+1) = hm;  lab{end+1} = sprintf('Feasible edge, m = %g kg', m_list(im));
    end
end
[Cw, hw] = contour(KK, CC, WORST, [1 1.1 1.25 1.5], ':', 'LineColor', [0.4 0.4 0.4]);
clabel(Cw, hw, 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
h(end+1) = hw;  lab{end+1} = 'Worst ratio to margin-adjusted limits';
h(end+1) = plot(k_best, c_best, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 13);
lab{end+1} = 'Selected design';

set(gca, 'XScale', axis_scale, 'YScale', axis_scale, 'Layer', 'top'); box on
xlim([min(k_vec) max(k_vec)]); ylim([min(c_vec) max(c_vec)]);
xlabel('Spring constant k [N/m]'); ylabel('Damper constant c [N\cdot{}s/m]');
title({sprintf('All masses: %d meet limits, %d meet margins', ...
      nnz(PASS_ALL), nnz(PASS_M_ALL)), ['Selected: ' sel_txt]});
legend(h, lab, 'Location', 'northeast');

%% ======================= PLOT: TIME HISTORY OF SELECTED DESIGN =======================
figure('Name', 'Selected design', 'Color', 'w');
for im = 1:nm
    [~, ~, ~, ~, z, am, DR] = run_case(m_list(im), k_best, c_best, u, a0, t, post, dt, g, ...
                                       stroke_def, accel_def, BR);
    subplot(3,1,1); hold on; plot(t, am/g, 'Color', cols(im,:), 'LineWidth', 1.3);
    subplot(3,1,2); hold on; plot(t, -z,   'Color', cols(im,:), 'LineWidth', 1.3);
    subplot(3,1,3); hold on; plot(t, DR,   'Color', cols(im,:), 'LineWidth', 1.3);
end
subplot(3,1,3); yline(BR.lim(1), 'r--', sprintf('+DR limit %.1f', BR.lim(1)));
yline(BR.lim(2), 'r--', sprintf('-DR limit %.1f', BR.lim(2)));
if beta_margin > 0
    yline(b_lim*BR.lim(1), 'r:');  yline(b_lim*BR.lim(2), 'r:');
end
grid on; ylabel(sprintf('DR%s (%s)', upper(brinkley_axis), crew_condition)); xlabel('t [s]');
subplot(3,1,1); plot(t, a/g, 'k:'); if isfinite(accel_max), yline(accel_max, 'r--'); end
if accel_margin > 0, yline(a_lim, 'r:'); end; grid on
ylabel('Accel reading [g]');
title(sprintf('k = %.0f N/m, c = %.0f N\\cdot{}s/m', k_best, c_best));
legend([arrayfun(@(m) sprintf('m = %g kg', m), m_list, 'UniformOutput', false), ...
        {'Capsule input'}], 'Location', 'northeast');
subplot(3,1,2); yline(stroke_max, 'r--', 'limit');
if stroke_margin > 0, yline(s_lim, 'r:', 'margin'); end; grid on
ylabel('Compression -z [m]');

%% ======================= LOCAL FUNCTIONS =======================
function [stroke, apk, imp, beta, z, am, DR] = run_case(m, k, c, u, a0, t, post, dt, g, ...
                                                       stroke_def, accel_def, BR)
% Exact discretization of the SDOF model, applied with filter().
% Simulates the deviation from equilibrium (zero initial state), then adds
% back the static deflection.
    Ac = [0 1; -k/m -c/m];
    B  = [0; -1];
    Ad = expm(Ac*dt);
    Bd = Ac \ ((Ad - eye(2))*B);
    tr = trace(Ad);
    den = [1, -tr, det(Ad)];
    M2  = Ad - tr*eye(2);

    Cz = [1 0];                 % output: z
    Ca = [-k/m -c/m];           % output: accelerometer reading on mass
    zdev = filter([0, Cz*Bd, Cz*M2*Bd], den, u);
    adev = filter([0, Ca*Bd, Ca*M2*Bd], den, u);

    z_static = -m*a0/k;
    z  = z_static + zdev;
    am = a0 + adev;

    switch stroke_def
        case 'total',   stroke = max(-z);
        case 'dynamic', stroke = max(-zdev);
        otherwise, error('stroke_def must be ''total'' or ''dynamic''');
    end
    switch accel_def
        case 'reading',   apk = max(abs(am))/g;
        case 'kinematic', apk = max(abs(am - a0))/g;
        otherwise, error('accel_def must be ''reading'' or ''kinematic''');
    end
    imp = trapz(t(post), abs(am(post) - a0))/g;

    % Brinkley: the seat (stroking mass) acceleration drives the occupant model
    if BR.on
        [beta, DR] = brinkley_dr(am, a0, dt, g, BR);
    else
        beta = 0;  DR = zeros(size(t));
    end
end

function [beta, DR] = brinkley_dr(am, a0, dt, g, BR)
% Brinkley Dynamic Response for one axis.
%   d'' + 2*zeta*wn*d' + wn^2*d = A(t),   DR(t) = wn^2*d(t)/g
% A(t) is the seat acceleration along the chosen axis. The model starts in
% equilibrium with the pre-impact input (DR = A(0)/g), so only the deviation
% from that state is filtered, exactly as in run_case.
    if BR.inc1g, A_in = BR.sgn*am;  A0 = BR.sgn*a0;
    else,        A_in = BR.sgn*(am - a0);  A0 = 0;
    end
    w  = BR.wn;  zt = BR.zeta;
    Ac = [0 1; -w^2 -2*zt*w];
    B  = [0; 1];
    Ad = expm(Ac*dt);
    Bd = Ac \ ((Ad - eye(2))*B);
    tr = trace(Ad);
    den = [1, -tr, det(Ad)];
    M2  = Ad - tr*eye(2);
    Cd  = [w^2 0];                         % output: wn^2 * d

    du = A_in - A0;
    du = [(du(1:end-1) + du(2:end))/2; du(end)];   % interval average
    DR = (A0 + filter([0, Cd*Bd, Cd*M2*Bd], den, du))/g;

    % Single axis: ratio to the limit that matches the sign of DR at each instant
    beta = max([max(DR)/BR.lim(1), min(DR)/BR.lim(2), 0]);
end

function add_zeta_lines(k_vec, c_vec, m, zetas)
% Lines of constant damping ratio: c = 2*zeta*sqrt(k*m)
    kk = linspace(min(k_vec), max(k_vec), 200);
    for zz = zetas
        cc = 2*zz*sqrt(kk*m);
        in = cc >= min(c_vec) & cc <= max(c_vec);
        if any(in)
            plot(kk(in), cc(in), '--', 'Color', [0.55 0.55 0.55], ...
                 'HandleVisibility', 'off');
            j = find(in, 1, 'last');
            text(kk(j), cc(j), sprintf(' \\zeta=%.1f', zz), 'FontSize', 8, ...
                 'Color', [0.45 0.45 0.45], 'VerticalAlignment', 'bottom', ...
                 'HorizontalAlignment', 'right');
        end
    end
end

