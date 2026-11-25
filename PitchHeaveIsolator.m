%% PitchHeaveIsolator.m
% Pitch-heave model of the isolated stretcher at splashdown, merged from
% SpaceAmbulancePitchHeaveModel.m (equal isolators), PitchHeave4Params.m
% (separate head/foot isolators, 4-D grid) and SpaceAmbulancePitchHeaveOptimizer.m
% (continuous polish). Revision 2026-11-20: each pulse shape runs at its own
% peak for the chosen pulse basis (isolator_envelope_peak.m; 'impulse' = 5 ms
% bounding pulse, floor mount). Revision 2026-11-19:
%   - inputs come from isolator_inputs.m (same file as the SDOF scripts)
%   - the pulse comes from isolator_make_input.m (same record as the SDOF scripts)
%   - Brinkley DR (X axis, supine) added at the PATIENT CG as two more states
%   - limits: stroke + Brinkley beta (legacy 10 g / 2 g*s only if use_legacy)
%   - self-check: with the CG mid-stretcher and equal isolators the model must
%     reduce to the SDOF model; compared with results_vibrations.mat if present
%
% Rigid stretcher + patient + bridge + equipment (mass m, pitch inertia J about
% the system CG) on two vertical spring-dampers, head end (z = 0) and foot end.
%
% Bond graph:  Sf: v_floor(t) at both isolator bases, Se: mg on 1:v_g,
%              I: m on 1:v_g, I: J on 1:w, TF: b and TF: a to the isolators,
%              C: 1/k_f, R: c_f at the head;  C: 1/k_r, R: c_r at the foot.
%              Brinkley body: unit mass driven (MSf, signal only) by the motion of
%              the stretcher point under the patient CG.
%
% Geometry (z from the head end):  b = z_cg - z_head,  a = z_foot - z_cg,
%   e = z_pat - z_cg (patient CG relative to the system CG)
%
% States  x = [p_m; p_J; q_f; q_r; d; d_dot]
%   v_f = v_g - b*w,  v_r = v_g + a*w
%   F_f = k_f q_f + c_f (v_floor - v_f),   F_r = k_r q_r + c_r (v_floor - v_r)
%   dp_m = F_f + F_r - m g,   dp_J = a F_r - b F_f
%   dq_f = v_floor - v_f,     dq_r = v_floor - v_r
%   felt acceleration at point z:  a(z) = (F_f + F_r)/m + (z - z_cg)*(a F_r - b F_f)/J
%   Brinkley:  d'' + 2 zeta wn d' + wn^2 d = a(z_pat) - g + gB,   DR = sgn*wn^2*d/g
%
% Base MATLAB only. Needs R2018b+ (local functions in a script).

clear; close all; clc;

%% ============================ INPUTS ==================================
IN = isolator_inputs();
g  = IN.g;
nP = numel(IN.m_list);
for ip = 1:nP
    pat(ip).m    = IN.m_list(ip);      %#ok<SAGROW>
    pat(ip).J    = IN.ph_J(ip);
    pat(ip).zcg  = IN.ph_zcg(ip);
    pat(ip).zpat = IN.ph_zpat(ip);
    pat(ip).name = sprintf('%.1f kg patient (M = %.1f kg)', IN.m_patient(ip), IN.m_list(ip));
end
G.zF = IN.ph_z_head;  G.zR = IN.ph_z_foot;  G.g = g;

% Brinkley
bt = IN.BTAB.(IN.brinkley_axis);
BR.on = IN.use_brinkley;  BR.wn = bt{1};  BR.zeta = bt{2};  BR.sgn = IN.brinkley_sign;
switch IN.crew_condition
    case 'healthy', BR.lim = bt{3};
    case 'ill',     BR.lim = bt{4};
end

% Limits
LIM.stroke = (1 - IN.stroke_margin)*IN.stroke_max;
LIM.beta   = (1 - IN.beta_margin);
LIM.legacy = IN.use_legacy;
LIM.acc    = (1 - IN.accel_margin)*IN.accel_max;
LIM.imp    = (1 - IN.impulse_margin)*IN.impulse_max;

% Inputs: one floor-velocity record per pulse shape (same generator as the SDOF scripts)
INr = IN;  INr.t_post = IN.ph_t_post;
for ii = 1:numel(IN.ph_pulses)
    INr.pulse = IN.ph_pulses{ii};
    INr.a_in_pk = isolator_envelope_peak(IN, INr.pulse);   % each shape at its own peak for the basis
    [t, acc, a0, t_imp, dt, T] = isolator_make_input(INr);
    dv = cumtrapz(t, acc - a0);
    PUL(ii).name = IN.ph_pulses{ii};                       %#ok<SAGROW>
    PUL(ii).t = t;  PUL(ii).dt = dt;  PUL(ii).a = acc;  PUL(ii).a0 = a0;
    PUL(ii).u = dv - dv(end);                              % floor velocity: -Vi before, 0 after
    PUL(ii).post = t >= t_imp;
    PUL(ii).apk = INr.a_in_pk;  PUL(ii).T = T;
end
BR.gB = g*IN.brinkley_1g + (g - PUL(1).a0)*(~IN.brinkley_1g);

fprintf('=== Pitch-heave inputs ===\n');
fprintf('Isolators at z = %.3f (head) and %.3f m (foot); Vi = %.2f m/s\n', G.zF, G.zR, IN.Vi);
for ii = 1:numel(PUL)
    fprintf('  pulse %-8s peak %.2f g, T = %.1f ms (%s basis)\n', PUL(ii).name, PUL(ii).apk, 1e3*PUL(ii).T, IN.pulse_basis);
end
for ip = 1:nP
    [b, a] = levers(pat(ip), G);
    fprintf('  %-32s J = %6.2f kg m^2, z_cg = %.4f, z_pat = %.4f m, b = %.4f, a = %.4f m\n', ...
            pat(ip).name, pat(ip).J, pat(ip).zcg, pat(ip).zpat, b, a);
end
fprintf('Pulses: %s;  limits: stroke %.2f m, beta %.2f (%s %s)%s\n\n', strjoin(IN.ph_pulses, ', '), ...
        LIM.stroke, LIM.beta, upper(IN.brinkley_axis), IN.crew_condition, ...
        ternary(LIM.legacy, sprintf(', %g g, %g g*s', LIM.acc, LIM.imp), ''));

%% ---------------- Self-check: reduction to the SDOF model ----------------
if exist('results_vibrations.mat', 'file')
    V = load('results_vibrations.mat');
    Lh = G.zR - G.zF;
    fprintf('Self-check vs SDOF (CG and patient mid-stretcher, equal isolators k/2, c/2):\n');
    INs = IN;  INs.pulse = 'halfsine';
    [ts, as, a0s, tis, dts] = isolator_make_input(INs);
    dvs = cumtrapz(ts, as - a0s);  us = dvs - dvs(end);
    for ip = 1:nP
        ps = pat(ip);  ps.zcg = G.zF + Lh/2;  ps.zpat = ps.zcg;
        D  = [V.chk.k/2; V.chk.k/2; V.chk.c/2; V.chk.c/2];
        H  = simulate_hist(ps, us, D, G, BR, dts);
        s1 = max(V.chk.stroke(:, ip));  s2 = max(max(H.qf), max(H.qr));
        a1 = max(abs(V.chk.am(:, ip)))/g;  a2 = max(abs(H.ag))/g;
        d1 = max(V.chk.DR(:, ip));  d2 = max(H.DR);
        fprintf('  m = %.2f kg: stroke %.4f vs %.4f m, peak %.3f vs %.3f g, DR peak %.3f vs %.3f (max rel. diff %.1e)\n', ...
                ps.m, s1, s2, a1, a2, d1, d2, max(abs([s2/s1, a2/a1, d2/d1] - 1)));
    end
    fprintf('\n');
end

%% ---------------- Equal isolators: 2-D grid ----------------
[KE, CE] = ndgrid(IN.ph_k_equal, IN.ph_c_equal);
De = [KE(:) KE(:) CE(:) CE(:)].';
fprintf('Equal isolators: %d designs\n', size(De, 2));
tic;  [We, Me] = eval_designs(De, pat, PUL, G, BR, LIM);  fprintf('  %.1f s\n', toc);
[wE, iE] = max(We);  bestEq = De(:, iE);
fprintf('  feasible %d / %d; best k = %.0f, c = %.0f per isolator, margin %.3f\n\n', ...
        nnz(We >= 0), numel(We), bestEq(1), bestEq(3), wE);

%% ---------------- Separate isolators: coarse 4-D grid around the equal optimum ----------------
kv = refine_axis(bestEq(1), IN.ph_sep_span(1)*bestEq(1), IN.ph_sep_n(1));
cv = refine_axis(bestEq(3), IN.ph_sep_span(2)*bestEq(3), IN.ph_sep_n(2));
[KF, KR, CF, CR] = ndgrid(kv, kv, cv, cv);
Ds = [KF(:) KR(:) CF(:) CR(:)].';
fprintf('Separate isolators: %d designs (k %.0f-%.0f, c %.0f-%.0f per isolator)\n', ...
        size(Ds, 2), kv(1), kv(end), cv(1), cv(end));
tic;  [Ws, ~] = eval_designs(Ds, pat, PUL, G, BR, LIM);  fprintf('  %.1f s\n', toc);
[wS, iS] = max(Ws);  bestSep = Ds(:, iS);

% Continuous polish (Nelder-Mead on log scale, base MATLAB)
if IN.ph_polish
    f = @(z) -eval_designs(exp(z(:)), pat, PUL, G, BR, LIM);
    zb = fminsearch(f, log(bestSep), optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
                                             'MaxFunEvals', 400));
    if -f(zb) > wS, bestSep = exp(zb(:));  wS = -f(zb); end
end
[~, Ms] = eval_designs(bestSep, pat, PUL, G, BR, LIM);
[~, Me1] = eval_designs(bestEq, pat, PUL, G, BR, LIM);

fprintf('\n=== Results ===\n');
fprintf('Best equal isolators:    k  = %6.0f N/m,             c  = %5.0f N*s/m        margin %+.2f%%\n', ...
        bestEq(1), bestEq(3), 100*wE);
fprintf('Best separate isolators: kf = %6.0f, kr = %6.0f N/m, cf = %5.0f, cr = %5.0f   margin %+.2f%%\n', ...
        bestSep, 100*wS);
fprintf('  head/foot ratios: k_f/k_r = %.3f, c_f/c_r = %.3f;  gain from separate tuning: %+.2f points\n', ...
        bestSep(1)/bestSep(2), bestSep(3)/bestSep(4), 100*(wS - wE));
print_table('Equal isolators', Me1, pat, PUL, bestEq, G);
print_table('Separate isolators', Ms, pat, PUL, bestSep, G);

% Independent check of the pitch coupling: Newton-Euler + ode45 (no bond graph,
% no discretization), heaviest patient, first pulse, best equal design
if IN.ph_ode45_check
    p = pat(end);  P1 = PUL(1);  D = bestEq;  [b, a] = levers(p, G);
    H = simulate_hist(p, P1.u, D, G, BR, P1.dt);
    x0 = equilibrium(p, D, G, BR, P1.u(1));
    zF = cumtrapz(P1.t, P1.u);
    f = @(tt, y) newton_euler(tt, y, P1.t, P1.u, zF, p.m, p.J, b, a, D, x0(3), x0(4), g);
    [~, Y] = ode45(f, P1.t, [0; 0; P1.u(1); 0], odeset('RelTol', 1e-9, 'AbsTol', 1e-12, 'MaxStep', P1.dt));
    qf = x0(3) + zF - (Y(:, 1) - b*Y(:, 2));  qr = x0(4) + zF - (Y(:, 1) + a*Y(:, 2));
    fprintf('\node45 Newton-Euler check (%s, %s): max |q_f| diff %.3g mm, max |q_r| diff %.3g mm (peaks %.1f / %.1f mm)\n', ...
            p.name, P1.name, 1e3*max(abs(H.qf(:) - qf)), 1e3*max(abs(H.qr(:) - qr)), 1e3*max(qf), 1e3*max(qr));
end

save('results_pitchheave.mat', 'IN', 'pat', 'bestEq', 'wE', 'bestSep', 'wS', 'Me1', 'Ms', 'De', 'We');

%% ---------------- Figures ----------------
% (1) equal-isolator design space: worst-case margin
figure('Name', 'Equal isolators: design space', 'Color', 'w');
imagesc(IN.ph_c_equal, IN.ph_k_equal, 100*reshape(We, size(KE)));  axis xy;  hold on
contour(IN.ph_c_equal, IN.ph_k_equal, reshape(We, size(KE)), [0 0], 'k-', 'LineWidth', 1.5);
plot(bestEq(3), bestEq(1), 'kp', 'MarkerSize', 13, 'MarkerFaceColor', 'w');
plot(bestSep(3:4), bestSep(1:2), 'ko', 'MarkerFaceColor', [0.6 0.6 0.6]);
colorbar;  caxis([-50 50]);
xlabel('c per isolator (N\cdots/m)');  ylabel('k per isolator (N/m)');
title({'Worst-case margin (%) over patients and pulses; black line = 0', ...
       'star: best equal; circles: best separate (head, foot)'});

% (2) time histories of the best equal design, one figure per pulse
cs = {[0.12 0.47 0.71], [0.84 0.15 0.16]};
for ii = 1:numel(PUL)
    figure('Name', ['Best equal design, ' PUL(ii).name], 'Color', 'w');
    tm = (PUL(ii).t - PUL(ii).t(find(PUL(ii).post, 1)))*1e3;
    for ip = 1:nP
        H = simulate_hist(pat(ip), PUL(ii).u, bestEq, G, BR, PUL(ii).dt);
        subplot(4,1,1); hold on; plot(tm, H.ap/g, 'Color', cs{ip}, 'DisplayName', pat(ip).name);
        subplot(4,1,2); hold on; plot(tm, H.qf, '-', 'Color', cs{ip}); plot(tm, H.qr, '--', 'Color', cs{ip});
        subplot(4,1,3); hold on; plot(tm, H.th*180/pi, 'Color', cs{ip});
        subplot(4,1,4); hold on; plot(tm, H.DR, 'Color', cs{ip});
    end
    subplot(4,1,1); grid on; ylabel('felt at patient CG (g)'); legend('Location', 'northeast');
    title(sprintf('%s pulse, k = %.0f N/m, c = %.0f N\\cdots/m per isolator', PUL(ii).name, bestEq(1), bestEq(3)));
    subplot(4,1,2); grid on; ylabel('compression (m)'); yline(LIM.stroke, 'r--');
    title('solid: head isolator, dashed: foot isolator');
    subplot(4,1,3); grid on; ylabel('pitch (deg, + foot up)');
    subplot(4,1,4); grid on; ylabel(sprintf('DR_%s', upper(IN.brinkley_axis))); xlabel('time from impact (ms)');
    yline(LIM.beta*BR.lim(1), 'r--');  yline(LIM.beta*BR.lim(2), 'r--');
end

%% ======================= local functions ==============================
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function v = refine_axis(center, halfwidth, n)
v = linspace(center - halfwidth, center + halfwidth, n);
v = v(v > 0);
end

function [b, a] = levers(p, G)
b = p.zcg - G.zF;  a = G.zR - p.zcg;
end

function [A, Bu, Bg, Cp, Dp] = ss6(p, G, BR, kf, kr, cf, cr)
% Continuous 6-state model. Inputs: floor velocity u, and the constant g.
[b, a] = levers(p, G);  m = p.m;  J = p.J;  e = p.zpat - p.zcg;
Ff = [-cf/m,  cf*b/J, kf, 0, 0, 0];   Ffu = cf;      % F_f = Ff*x + Ffu*u
Fr = [-cr/m, -cr*a/J, 0, kr, 0, 0];   Fru = cr;
Cp = (Ff + Fr)/m + e*(a*Fr - b*Ff)/J;                % felt accel at the patient CG
Dp = (Ffu + Fru)/m + e*(a*Fru - b*Ffu)/J;
w = BR.wn;  z = BR.zeta;
A  = [Ff + Fr;
      a*Fr - b*Ff;
      -1/m,  b/J, 0, 0, 0, 0;
      -1/m, -a/J, 0, 0, 0, 0;
       0,    0,   0, 0, 0, 1;
       Cp + [0 0 0 0 -w^2 -2*z*w]];
Bu = [Ffu + Fru; a*Fru - b*Ffu; 1; 1; 0; Dp];
Bg = [-m; 0; 0; 0; 0; (BR.gB - G.g)/G.g];           % times the input g
end

function x0 = equilibrium(p, D, G, BR, u0)
[b, a] = levers(p, G);  n = size(D, 2);
x0 = [p.m*u0*ones(1, n); zeros(1, n); ...
      p.m*G.g*a./(D(1, :)*(a + b)); p.m*G.g*b./(D(2, :)*(a + b)); ...
      BR.gB/BR.wn^2*ones(1, n); zeros(1, n)];
end

function [Ad, Bud, Bgd] = zoh(A, Bu, Bg, dt)
M = zeros(8);  M(1:6, 1:6) = A;  M(1:6, 7) = Bu;  M(1:6, 8) = Bg;
E = expm(M*dt);
Ad = E(1:6, 1:6);  Bud = E(1:6, 7);  Bgd = E(1:6, 8);
end

function [W, M] = eval_designs(D, pat, PUL, G, BR, LIM)
% Worst-case margin over patients and pulses for each design (column of D)
n = size(D, 2);  W = inf(1, n);
M = struct('stroke', {}, 'ext', {}, 'beta', {}, 'peak', {}, 'imp', {}, 'pitch', {});
for ip = 1:numel(pat)
    for ii = 1:numel(PUL)
        R = sweep_metrics(pat(ip), PUL(ii), D, G, BR);
        M(ip, ii) = R;
        r = [R.stroke/LIM.stroke; BR.on*R.beta/LIM.beta];
        if LIM.legacy, r = [r; R.peak/LIM.acc; R.imp/LIM.imp]; end %#ok<AGROW>
        W = min(W, 1 - max(r, [], 1));
    end
end
end

function R = sweep_metrics(p, P, D, G, BR)
% Vectorized over designs; metrics accumulated on the fly
n = size(D, 2);  dt = P.dt;  [b, a] = levers(p, G);
Adf = zeros(36, n);  Bu = zeros(6, n);  Bg = zeros(6, n);
for i = 1:n
    [A, Bu1, Bg1] = ss6(p, G, BR, D(1, i), D(2, i), D(3, i), D(4, i));
    [Ad, Bu(:, i), Bg(:, i)] = zoh(A, Bu1, Bg1, dt);
    Adf(:, i) = Ad(:);
end
x  = equilibrium(p, D, G, BR, P.u(1));
ub = [(P.u(1:end-1) + P.u(2:end))/2; P.u(end)];          % interval-average velocity (2nd order)
kf = D(1, :);  kr = D(2, :);  cf = D(3, :);  cr = D(4, :);
R.stroke = -inf(1, n);  R.ext = inf(1, n);  drmx = -inf(1, n);  drmn = inf(1, n);
R.peak = zeros(1, n);  imp = zeros(3, n);  R.pitch = zeros(1, n);
for j = 1:numel(P.u)
    if j > 1
        xn = Bu*ub(j - 1) + Bg*G.g;
        for c = 1:6, xn = xn + Adf(6*c-5:6*c, :).*x(c, :); end
        x = xn;
    end
    vg = x(1, :)/p.m;  w = x(2, :)/p.J;
    Ff = kf.*x(3, :) + cf.*(P.u(j) - (vg - b*w));
    Fr = kr.*x(4, :) + cr.*(P.u(j) - (vg + a*w));
    ag = (Ff + Fr)/p.m;  al = (a*Fr - b*Ff)/p.J;
    accs = [ag; ag - b*al; ag + a*al]/G.g;                % felt: CG, head end, foot end
    R.peak = max(R.peak, max(abs(accs), [], 1));
    if P.post(j), imp = imp + abs(accs - 1)*dt; end
    R.stroke = max(R.stroke, max(x(3, :), x(4, :)));
    R.ext = min(R.ext, min(x(3, :), x(4, :)));
    R.pitch = max(R.pitch, abs(x(3, :) - x(4, :))/(a + b));
    DR = BR.sgn*BR.wn^2*x(5, :)/G.g;
    drmx = max(drmx, DR);  drmn = min(drmn, DR);
end
R.imp  = max(imp, [], 1);
R.beta = max([drmx/BR.lim(1); drmn/BR.lim(2); zeros(1, n)], [], 1);
end

function H = simulate_hist(p, u, D, G, BR, dt)
[b, a] = levers(p, G);
[A, Bu, Bg] = ss6(p, G, BR, D(1), D(2), D(3), D(4));
[Ad, Bud, Bgd] = zoh(A, Bu, Bg, dt);
ub = [(u(1:end-1) + u(2:end))/2; u(end)];
X = zeros(6, numel(u));  X(:, 1) = equilibrium(p, D, G, BR, u(1));
for j = 2:numel(u), X(:, j) = Ad*X(:, j - 1) + Bud*ub(j - 1) + Bgd*G.g; end
vg = X(1, :)/p.m;  w = X(2, :)/p.J;  uu = u(:).';
Ff = D(1)*X(3, :) + D(3)*(uu - (vg - b*w));
Fr = D(2)*X(4, :) + D(4)*(uu - (vg + a*w));
H.ag = (Ff + Fr)/p.m;  al = (a*Fr - b*Ff)/p.J;
H.ap = H.ag + (p.zpat - p.zcg)*al;
H.qf = X(3, :);  H.qr = X(4, :);  H.th = (X(3, :) - X(4, :))/(a + b);
H.DR = BR.sgn*BR.wn^2*X(5, :)/G.g;
end

function dy = newton_euler(tt, y, t, u, zFloor, m, J, b, a, D, qf0, qr0, g)
% y = [z; theta; z_dot; theta_dot] measured from the initial equilibrium
vb = interp1(t, u, tt);  zb = interp1(t, zFloor, tt);
Ff = D(1)*(qf0 + zb - (y(1) - b*y(2))) + D(3)*(vb - (y(3) - b*y(4)));
Fr = D(2)*(qr0 + zb - (y(1) + a*y(2))) + D(4)*(vb - (y(3) + a*y(4)));
dy = [y(3); y(4); (Ff + Fr)/m - g; (a*Fr - b*Ff)/J];
end

function print_table(ttl, M, pat, PUL, D, G)
fprintf('\n%s (k_f %.0f, k_r %.0f N/m; c_f %.0f, c_r %.0f N*s/m):\n', ttl, D);
for ip = 1:numel(pat)
    [b, a] = levers(pat(ip), G);
    fprintf('  %-32s static sag head %.3f m, foot %.3f m\n', pat(ip).name, ...
            pat(ip).m*G.g*a/(D(1)*(a + b)), pat(ip).m*G.g*b/(D(2)*(a + b)));
    for ii = 1:numel(PUL)
        R = M(ip, ii);
        fprintf('    %-9s stroke %.3f m, beta %.3f, peak felt %5.2f g, pitch %.2f deg, min compression %.3f m\n', ...
                PUL(ii).name, R.stroke, R.beta, R.peak, R.pitch*180/pi, R.ext);
    end
end
end
