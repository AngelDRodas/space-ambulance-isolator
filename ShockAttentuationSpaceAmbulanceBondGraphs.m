%% ShockAttentuationSpaceAmbulanceBondGraphs.m  (was splashdown_sweep_bondgraph.m)
% Revision 2026-11-19: inputs moved to isolator_inputs.m (shared with the
% vibrations script); pulse generation moved to isolator_make_input.m;
% Brinkley DR added as two more bond-graph states (see below); saves
% results_bondgraph.mat for compare_models.m.
%
% Bond-graph formulation of the splashdown isolated-stretcher sweep.
% Same study, inputs, metrics, selection rule and plots as the vibrations
% script (ShockAttentuationSpaceAmbulanceBrinkley.m), but the dynamics come
% from the bond graph state equations and are integrated directly
% (fixed-step RK4), not from the vibrations EOM/filter().
% The two scripts are independent paths to the same answer, for cross-validation.
%
% Isolator bond graph (Sf: base velocity, Se: gravity):
%
%   Sf --1--> 1_vb --2--> 0 --6--> 1_vm --8--> I:m
%                         |          ^
%                         3          7
%                         |          |
%                         |          Se:-mg
%                         v
%            R:b <--4-- 1_(vb-vm) --5--> C:1/k
%
% States:  q5 = spring compression  (x_b - x_m, positive = compressed)
%          p8 = momentum of the mass (m*v_m, positive UP)
%
%   q5_dot = f5 = v_b - p8/m
%   p8_dot = e8 = k*q5 + b*(v_b - p8/m) - m*g
%
% Brinkley body sub-graph (driven by the stretcher, no load back on it):
%
%   MSf(v_m) --> 1_vm' --> 0 --> 1_vB --> I:1        (unit body mass)
%                          |       ^
%                          |       Se:-gB
%                          v
%      R:2*zeta*wn <-- 1_(vm-vB) --> C:1/wn^2
%
% The Brinkley model is a unit-mass spring-damper sitting on the stretcher.
% It is driven by a MODULATED flow source equal to the stretcher velocity
% v_m = p8/m (a signal, no power returned), because Brinkley's body mass is
% already part of m and must not load the isolator a second time.
% States:  qB = body "spinal/chest" compression (positive = compressed)
%          pB = body momentum per unit mass (= v_B)
%   qB_dot = v_m - pB
%   pB_dot = wn^2*qB + 2*zeta*wn*(v_m - pB) - gB
%   DR(t)  = sgn * wn^2 * qB / g
% gB = g keeps the 1 g in the input (DR = 1 at rest, brinkley_1g = true);
% gB = g - a0 removes it (DR = 0 at rest). Eliminating vB gives exactly
%   qB'' + 2*zeta*wn*qB' + wn^2*qB = (stretcher reading) - g + gB,
% the same Brinkley equation the vibrations script filters.
%
% Input: the same capsule accelerometer reading a(t) (specific force) used by
% the reference script. The flow source is v_b(t), so v_b is carried as a
% state with v_b_dot = a(t) - g; gravity enters separately through Se.
%   Accelerometer reading on the mass = p8_dot/m + g = (e5 + e4)/m
%   Static equilibrium: k*q5 = m*a0  (same as z(0) = -m*a0/k, z = -q5)
%
% Metrics match the reference exactly:
%   stroke  'total' = max q5,  'dynamic' = max(q5 - m*a0/k)
%   accel   'reading' = max|reading|, 'kinematic' = max|reading - a0|
%   impulse integral of |reading - a0| dt from impact on [g*s]
%   beta    max DR / limit (limit picked by sign)
%
% All (k,b) designs for one mass are integrated together as vectors, so the
% sweep needs no per-design loop. No toolboxes required.

clear; clc; close all

%% ======================= INPUTS (shared) =======================
% All inputs live in isolator_inputs.m. Edit that file, not this one.
IN = isolator_inputs();   % (named IN: P is used below for pass masks)
g = IN.g;
Vi = IN.Vi;  pulse = IN.pulse;  a_in_pk = IN.a_in_pk;
stroke_max = IN.stroke_max;  stroke_def = IN.stroke_def;  accel_def = IN.accel_def;
if IN.use_legacy
    accel_max = IN.accel_max;  impulse_max = IN.impulse_max;
else
    accel_max = Inf;           impulse_max = Inf;    % legacy limits off
end
stroke_margin = IN.stroke_margin;  accel_margin = IN.accel_margin;
impulse_margin = IN.impulse_margin;  beta_margin = IN.beta_margin;
selection = IN.selection;
m_list = IN.m_list;  k_vec = IN.k_vec;  b_vec = IN.c_vec;   % R:b = damper c
axis_scale = IN.axis_scale;  zeta_guides = IN.zeta_guides;
out_all    = 'sweep_results_bondgraph.csv';
out_common = 'feasible_all_masses_bondgraph.csv';
fprintf('Isolated mass: %.2f kg package + patient %s kg -> m = %s kg\n', ...
        IN.m_package, mat2str(IN.m_patient), mat2str(m_list));

%% ======================= INPUT PROFILE =======================
[t, a, a0, t_imp, dt, T, post] = isolator_make_input(IN);
vb0  = -trapz(t, a - g);        % initial base velocity: capsule ends at rest
% (Only relative motion matters, so vb0 doesn't change any metric. It just
%  makes v_b and v_m physically meaningful in the time-history plot.)

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
bt = IN.BTAB.(IN.brinkley_axis);
BR.on   = IN.use_brinkley;
BR.wn   = bt{1};  BR.zeta = bt{2};
switch IN.crew_condition
    case 'healthy', BR.lim = bt{3};
    case 'ill',     BR.lim = bt{4};
    otherwise, error('crew_condition must be ''healthy'' or ''ill''');
end
BR.sgn = IN.brinkley_sign;
if IN.brinkley_1g, BR.gB = g; else, BR.gB = g - a0; end   % Se on the body mass
fprintf('Brinkley %s-axis, %s crew: wn = %.1f rad/s, zeta = %.3f, DR limits +%.1f / %.1f\n\n', ...
        upper(IN.brinkley_axis), IN.crew_condition, BR.wn, BR.zeta, BR.lim(1), BR.lim(2));

%% ======================= SWEEP =======================
nm = numel(m_list);  nk = numel(k_vec);  nb = numel(b_vec);
[KK, BB] = meshgrid(k_vec, b_vec);        % rows: b, columns: k

STROKE = zeros(nb, nk, nm);
ACC    = zeros(nb, nk, nm);
IMP    = zeros(nb, nk, nm);
BETA   = zeros(nb, nk, nm);

tic
for im = 1:nm
    [s, ac, ip, be] = sim_bg(m_list(im), KK(:), BB(:), a, a0, t, post, g, ...
                             stroke_def, accel_def, vb0, BR);
    STROKE(:,:,im) = reshape(s,  nb, nk);
    ACC(:,:,im)    = reshape(ac, nb, nk);
    IMP(:,:,im)    = reshape(ip, nb, nk);
    BETA(:,:,im)   = reshape(be, nb, nk);
end
fprintf('Swept %d designs x %d masses in %.1f s (bond graph, RK4)\n\n', nk*nb, nm, toc);
if ~BR.on, BETA(:) = 0; end

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

RATIO = max(cat(4, STROKE/s_lim, ACC/a_lim, IMP/i_lim, BETA/b_lim), [], 4);
WORST = max(RATIO, [], 3);

%% ======================= SUMMARY =======================
for im = 1:nm
    P = PASS(:,:,im);
    fprintf('m = %g kg: %d of %d designs meet limits, %d also meet margins\n', ...
            m_list(im), nnz(P), numel(P), nnz(PASS_M(:,:,im)));
    if any(P(:))
        fprintf('   k range: %.0f - %.0f N/m\n',   min(KK(P)), max(KK(P)));
        fprintf('   b range: %.0f - %.0f N*s/m\n', min(BB(P)), max(BB(P)));
    end
end

S_W   = max(STROKE, [], 3);
OTHER = max(max(max(ACC/a_lim, IMP/i_lim), BETA/b_lim), [], 3);
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
[~, ibest] = min(obj(:));
k_best = KK(ibest);  b_best = BB(ibest);
fprintf('\nCommon to all masses: %d meet limits, %d also meet margins\n', ...
        nnz(PASS_ALL), nnz(PASS_M_ALL));
fprintf('Selected design (%s):\n  k = %.0f N/m, b = %.0f N*s/m\n', ...
        sel_txt, k_best, b_best);
fprintf('Margin = remaining fraction of each hard limit (negative = exceeded)\n');
fprintf('%8s %7s %10s %8s %9s %8s %12s %8s %7s %8s %8s\n', 'm [kg]', 'zeta', ...
        'stroke[m]', 'margin', 'accel[g]', 'margin', 'impulse[g*s]', 'margin', ...
        'beta', 'margin', 'DR peak');
for im = 1:nm
    [s, ac, ip, be, H] = sim_bg(m_list(im), k_best, b_best, a, a0, t, post, g, ...
                                stroke_def, accel_def, vb0, BR);
    fprintf('%8g %7.3f %10.3f %7.0f%% %9.2f %7.0f%% %12.3f %7.0f%% %7.3f %7.0f%% %8.2f\n', m_list(im), ...
            b_best/(2*sqrt(k_best*m_list(im))), ...
            s,  100*(1 - s/stroke_max), ac, 100*(1 - ac/accel_max), ...
            ip, 100*(1 - ip/impulse_max), be, 100*(1 - be), max(H.DR));
end

%% ======================= CSV OUTPUT =======================
% Column names match the reference CSV (c_Nspm = b).
MM3 = repmat(reshape(m_list, 1, 1, nm), nb, nk);
KK3 = repmat(KK, 1, 1, nm);
BB3 = repmat(BB, 1, 1, nm);
ZZ3 = BB3 ./ (2*sqrt(KK3 .* MM3));
results = table(MM3(:), KK3(:), BB3(:), ZZ3(:), STROKE(:), ACC(:), IMP(:), BETA(:), ...
    PASS(:), PASS_M(:), RATIO(:), ...
    'VariableNames', {'mass_kg','k_Npm','c_Nspm','zeta','stroke_m', ...
                      'accel_g','impulse_gs','brinkley_beta','pass_limits','pass_margins','margin_ratio'});
writetable(results, out_all);

if any(PASS_ALL(:))
    common = table(KK(PASS_ALL), BB(PASS_ALL), PASS_M_ALL(PASS_ALL), WORST(PASS_ALL), ...
        'VariableNames', {'k_Npm','c_Nspm','meets_margins','worst_margin_ratio'});
    common = sortrows(common, 'worst_margin_ratio');
    writetable(common, out_common);
end

% MAT file for compare_models.m: full sweep plus time histories of the
% check design (same k, c in both models).
[~, ~, ~, ~, Hc] = sim_bg_multi(m_list, IN.k_check, IN.c_check, a, a0, t, post, g, ...
                               stroke_def, accel_def, vb0, BR);
R = struct('model', 'bond graph (RK4)', 'IN', IN, 'm_list', m_list, ...
           'k_vec', k_vec, 'c_vec', b_vec, 'STROKE', STROKE, 'ACC', ACC, ...
           'IMP', IMP, 'BETA', BETA, 'PASS', PASS, 'k_best', k_best, ...
           'c_best', b_best, 't', t, 'a', a, 'dt', dt, ...
           'chk', struct('k', IN.k_check, 'c', IN.c_check, 'stroke', Hc.q5, ...
                         'am', Hc.am, 'DR', Hc.DR));
save('results_bondgraph.mat', '-struct', 'R');


%% ======================= PLOTS: EACH MASS =======================
for im = 1:nm
    m = m_list(im);
    figure('Name', sprintf('Bond graph: m = %g kg', m), 'Color', 'w'); hold on
    P  = PASS(:,:,im);
    PM = PASS_M(:,:,im);
    hP = plot([KK(P & ~PM); NaN], [BB(P & ~PM); NaN], 'o', 'MarkerSize', 4, ...
              'MarkerFaceColor', [0.65 0.85 0.65], 'MarkerEdgeColor', 'none');
    hM = plot([KK(PM); NaN], [BB(PM); NaN], 'o', 'MarkerSize', 4, ...
              'MarkerFaceColor', [0.1 0.5 0.1], 'MarkerEdgeColor', 'none');
    if stroke_margin > 0
        contour(KK, BB, STROKE(:,:,im), [s_lim s_lim], '--', ...
                'LineColor', [0.85 0.2 0.2], 'HandleVisibility', 'off');
    end
    if accel_margin > 0
        contour(KK, BB, ACC(:,:,im), [a_lim a_lim], '--', ...
                'LineColor', [0.15 0.35 0.8], 'HandleVisibility', 'off');
    end
    [~, hS] = contour(KK, BB, STROKE(:,:,im), [stroke_max stroke_max], ...
                      'LineColor', [0.85 0.2 0.2], 'LineWidth', 1.6);
    if isfinite(accel_max)
        [~, hA] = contour(KK, BB, ACC(:,:,im), [accel_max accel_max], ...
                          'LineColor', [0.15 0.35 0.8], 'LineWidth', 1.6);
        [~, hI] = contour(KK, BB, IMP(:,:,im), [impulse_max impulse_max], ...
                          'LineColor', [0.9 0.55 0.1], 'LineWidth', 1.6);
    else            % legacy limits off: placeholders keep the legend intact
        hA = plot(NaN, NaN, 'Color', [0.15 0.35 0.8], 'LineWidth', 1.6);
        hI = plot(NaN, NaN, 'Color', [0.9 0.55 0.1], 'LineWidth', 1.6);
    end
    if BR.on
        [~, hBr] = contour(KK, BB, BETA(:,:,im), [1 1], ...
                           'LineColor', [0.55 0.2 0.7], 'LineWidth', 1.6);
    else
        hBr = plot(NaN, NaN, 'Color', [0.55 0.2 0.7], 'LineWidth', 1.6);
    end
    add_zeta_lines(k_vec, b_vec, m, zeta_guides);

    [~, ibm] = min(reshape(RATIO(:,:,im), [], 1));
    hB = plot(KK(ibm), BB(ibm), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 13);

    set(gca, 'XScale', axis_scale, 'YScale', axis_scale, 'Layer', 'top'); box on
    xlim([min(k_vec) max(k_vec)]); ylim([min(b_vec) max(b_vec)]);
    xlabel('Spring constant k [N/m]'); ylabel('Damper constant b [N\cdot{}s/m]');
    title(sprintf('Bond graph, m = %g kg: %d meet limits, %d meet margins', ...
                  m, nnz(P), nnz(PM)));
    legend([hP hM hS hA hI hBr hB], {'Meets limits only', 'Meets limits + margins', ...
        sprintf('Stroke = %.2f m', stroke_max), ...
        sprintf('Accel = %g g', accel_max), ...
        sprintf('Impulse = %g g\\cdot{}s', impulse_max), ...
        sprintf('Brinkley \\beta = 1 (%s, %s)', upper(IN.brinkley_axis), IN.crew_condition), ...
        'Best for this mass'}, 'Location', 'northeast');
end

%% ======================= PLOT: ALL MASSES =======================
figure('Name', 'Bond graph: common solution space', 'Color', 'w'); hold on
cols = lines(nm);
h = gobjects(0);  lab = {};
LO = PASS_ALL & ~PASS_M_ALL;
h(end+1) = plot([KK(LO); NaN], [BB(LO); NaN], 'o', 'MarkerSize', 4, ...
                'MarkerFaceColor', [0.65 0.85 0.65], 'MarkerEdgeColor', 'none');
lab{end+1} = sprintf('All masses, limits only (%d)', nnz(LO));
h(end+1) = plot([KK(PASS_M_ALL); NaN], [BB(PASS_M_ALL); NaN], 'o', 'MarkerSize', 4, ...
                'MarkerFaceColor', [0.1 0.5 0.1], 'MarkerEdgeColor', 'none');
lab{end+1} = sprintf('All masses, limits + margins (%d)', nnz(PASS_M_ALL));
for im = 1:nm
    if any(any(PASS(:,:,im)))
        [~, hm] = contour(KK, BB, double(PASS(:,:,im)), [0.5 0.5], ...
                          'LineColor', cols(im,:), 'LineWidth', 1.8);
        h(end+1) = hm;  lab{end+1} = sprintf('Feasible edge, m = %g kg', m_list(im));
    end
end
[Cw, hw] = contour(KK, BB, WORST, [1 1.1 1.25 1.5], ':', 'LineColor', [0.4 0.4 0.4]);
clabel(Cw, hw, 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
h(end+1) = hw;  lab{end+1} = 'Worst ratio to margin-adjusted limits';
h(end+1) = plot(k_best, b_best, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 13);
lab{end+1} = 'Selected design';

set(gca, 'XScale', axis_scale, 'YScale', axis_scale, 'Layer', 'top'); box on
xlim([min(k_vec) max(k_vec)]); ylim([min(b_vec) max(b_vec)]);
xlabel('Spring constant k [N/m]'); ylabel('Damper constant b [N\cdot{}s/m]');
title({sprintf('Bond graph, all masses: %d meet limits, %d meet margins', ...
      nnz(PASS_ALL), nnz(PASS_M_ALL)), ['Selected: ' sel_txt]});
legend(h, lab, 'Location', 'northeast');

%% ======================= PLOT: TIME HISTORY OF SELECTED DESIGN =======================
figure('Name', 'Bond graph: selected design', 'Color', 'w');
for im = 1:nm
    m = m_list(im);
    [~, ~, ~, ~, H] = sim_bg(m, k_best, b_best, a, a0, t, post, g, ...
                             stroke_def, accel_def, vb0, BR);
    subplot(4,1,1); hold on; plot(t, H.am/g,  'Color', cols(im,:), 'LineWidth', 1.3);
    subplot(4,1,2); hold on; plot(t, H.q5,    'Color', cols(im,:), 'LineWidth', 1.3);
    subplot(4,1,3); hold on; plot(t, H.p8/m,  'Color', cols(im,:), 'LineWidth', 1.3);
    subplot(4,1,4); hold on; plot(t, H.DR,    'Color', cols(im,:), 'LineWidth', 1.3);
end
subplot(4,1,1); plot(t, a/g, 'k:'); if isfinite(accel_max), yline(accel_max, 'r--'); end
if accel_margin > 0, yline(a_lim, 'r:'); end; grid on
ylabel('Accel reading [g]');
title(sprintf('Bond graph: k = %.0f N/m, b = %.0f N\\cdot{}s/m', k_best, b_best));
legend([arrayfun(@(m) sprintf('m = %g kg', m), m_list, 'UniformOutput', false), ...
        {'Capsule input'}], 'Location', 'northeast');
subplot(4,1,2); yline(stroke_max, 'r--', 'limit');
if stroke_margin > 0, yline(s_lim, 'r:', 'margin'); end; grid on
ylabel('Compression q_5 [m]');
subplot(4,1,3); plot(t, H.vb, 'k:'); grid on
ylabel('Velocity [m/s]');
legend([arrayfun(@(m) sprintf('v_m = p_8/m, m = %g kg', m), m_list, ...
        'UniformOutput', false), {'v_b (capsule)'}], 'Location', 'southeast');
subplot(4,1,4); yline(BR.lim(1), 'r--', sprintf('+DR limit %.1f', BR.lim(1)));
yline(BR.lim(2), 'r--', sprintf('-DR limit %.1f', BR.lim(2)));
if beta_margin > 0, yline(b_lim*BR.lim(1), 'r:'); yline(b_lim*BR.lim(2), 'r:'); end
grid on; ylabel(sprintf('DR%s (%s)', upper(IN.brinkley_axis), IN.crew_condition)); xlabel('t [s]');

%% ======================= LOCAL FUNCTIONS =======================
function [stroke, apk, imp, beta, H] = sim_bg(m, K, B, a, a0, t, post, g, ...
                                              stroke_def, accel_def, vb0, BR)
% Integrates the bond graph state equations for every (K,B) pair at once.
% K, B: vectors of spring and damper constants (one entry per design).
% m may be a scalar or a vector matching K (one mass per design).
% Fixed-step RK4 on the sample grid; the input is linearly interpolated
% between samples (midpoint value used for the RK4 middle stages).
    K = K(:);  B = B(:);  m = m(:);
    nt = numel(t);  nd = numel(K);

    % Initial state: static equilibrium under the pre-impact reading a0
    q  = m*a0 ./ K;                   % q5:  k*q5 = m*a0
    vb = vb0;                         % base velocity (same for every design)
    p  = m.*vb.*ones(nd, 1);          % p8:  mass moves with the base
    % Brinkley body at rest on the stretcher: wn^2*qB = a0 - g + gB
    qB = (a0 - g + BR.gB)/BR.wn^2 * ones(nd, 1);
    pB = vb*ones(nd, 1);

    switch stroke_def
        case 'total',   s_off = 0;
        case 'dynamic', s_off = m*a0 ./ K;    % static compression
        otherwise, error('stroke_def must be ''total'' or ''dynamic''');
    end
    switch accel_def
        case 'reading',   a_off = 0;
        case 'kinematic', a_off = a0;
        otherwise, error('accel_def must be ''reading'' or ''kinematic''');
    end

    am   = reading(q, p, vb, m, K, B);
    DR   = BR.sgn*BR.wn^2*qB/g;
    smax = q - s_off;
    amax = abs(am - a_off);
    drmx = DR;  drmn = DR;
    imp  = zeros(nd, 1);
    prev = abs(am - a0);
    ip   = find(post, 1);

    keep = nargout > 4;
    if keep
        H.q5 = zeros(nt, nd);  H.p8 = zeros(nt, nd);
        H.vb = zeros(nt, 1);   H.am = zeros(nt, nd);  H.DR = zeros(nt, nd);
        H.q5(1,:) = q.';  H.p8(1,:) = p.';  H.vb(1) = vb;  H.am(1,:) = am.';
        H.DR(1,:) = DR.';
    end

    for n = 1:nt-1
        h  = t(n+1) - t(n);
        a1 = a(n);  a3 = a(n+1);  a2 = 0.5*(a1 + a3);
        X  = {q, p, vb, qB, pB};
        k1 = bg_rhs(X,                   a1, m, K, B, g, BR);
        k2 = bg_rhs(axpy(X, 0.5*h, k1),  a2, m, K, B, g, BR);
        k3 = bg_rhs(axpy(X, 0.5*h, k2),  a2, m, K, B, g, BR);
        k4 = bg_rhs(axpy(X, h, k3),      a3, m, K, B, g, BR);
        for j = 1:5
            X{j} = X{j} + (h/6)*(k1{j} + 2*k2{j} + 2*k3{j} + k4{j});
        end
        [q, p, vb, qB, pB] = X{:};

        am   = reading(q, p, vb, m, K, B);
        DR   = BR.sgn*BR.wn^2*qB/g;
        smax = max(smax, q - s_off);
        amax = max(amax, abs(am - a_off));
        drmx = max(drmx, DR);  drmn = min(drmn, DR);
        cur  = abs(am - a0);
        if n >= ip
            imp = imp + 0.5*h*(prev + cur);       % trapezoid, impact onward
        end
        prev = cur;

        if keep
            H.q5(n+1,:) = q.';  H.p8(n+1,:) = p.';
            H.vb(n+1)   = vb;   H.am(n+1,:) = am.';  H.DR(n+1,:) = DR.';
        end
    end
    stroke = smax;
    apk    = amax / g;
    imp    = imp  / g;
    % Single axis: ratio to the limit that matches the sign of DR
    beta   = max([drmx/BR.lim(1), drmn/BR.lim(2), zeros(nd,1)], [], 2);
end

function [stroke, apk, imp, beta, H] = sim_bg_multi(m_list, k, b, a, a0, t, post, g, ...
                                                    stroke_def, accel_def, vb0, BR)
% One design (k, b), every mass at once (used for the saved check design).
    nm = numel(m_list);
    [stroke, apk, imp, beta, H] = sim_bg(m_list(:), k*ones(nm,1), b*ones(nm,1), ...
                                         a, a0, t, post, g, stroke_def, accel_def, vb0, BR);
end

function Y = axpy(X, h, D)
% Y = X + h*D for the cell array of states
    Y = X;
    for j = 1:numel(X), Y{j} = X{j} + h*D{j}; end
end

function D = bg_rhs(X, ab, m, K, B, g, BR)
% Bond graph constitutive and junction relations (bond numbers as drawn).
    [q5, p8, vb, qB, pB] = X{:};
    % --- isolator ---
    f8  = p8 ./ m;           % I:m        flow from momentum
    f6  = f8;                % 1_vm       common flow
    f2  = vb;                % Sf, 1_vb   f1 = f2 = v_b
    f3  = f2 - f6;           % 0-junction flow balance
    f4  = f3;  f5 = f3;      % 1-junction common flow
    e5  = K .* q5;           % C:1/k
    e4  = B .* f4;           % R:b
    e3  = e4 + e5;           % 1-junction effort sum
    e6  = e3;                % 0-junction common effort
    e7  = -m*g;              % Se: gravity
    e8  = e6 + e7;           % 1_vm effort sum
    % --- Brinkley body (unit mass), MSf = v_m: signal only, no load back ---
    fB  = pB;                        % I:1   body velocity
    fr  = f8 - fB;                   % 0-junction: relative velocity v_m - v_B
    eC  = BR.wn^2 * qB;              % C:1/wn^2
    eR  = 2*BR.zeta*BR.wn * fr;      % R:2*zeta*wn
    eI  = eC + eR - BR.gB;           % 1_vB effort sum with Se:-gB
    D = {f5, e8, ab - g, fr, eI};    % q5', p8', vb', qB', pB'
end

function am = reading(q5, p8, vb, m, K, B)
% Accelerometer reading on the mass (specific force) = p8_dot/m + g = (e5+e4)/m
    am = (K.*q5 + B.*(vb - p8./m)) ./ m;
end

function add_zeta_lines(k_vec, b_vec, m, zetas)
% Lines of constant damping ratio: b = 2*zeta*sqrt(k*m)
    kk = linspace(min(k_vec), max(k_vec), 200);
    for zz = zetas
        bb = 2*zz*sqrt(kk*m);
        in = bb >= min(b_vec) & bb <= max(b_vec);
        if any(in)
            plot(kk(in), bb(in), '--', 'Color', [0.55 0.55 0.55], ...
                 'HandleVisibility', 'off');
            j = find(in, 1, 'last');
            text(kk(j), bb(j), sprintf(' \\zeta=%.1f', zz), 'FontSize', 8, ...
                 'Color', [0.45 0.45 0.45], 'VerticalAlignment', 'bottom', ...
                 'HorizontalAlignment', 'right');
        end
    end
end
