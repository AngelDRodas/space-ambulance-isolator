function [t, a, a0, t_imp, dt, T, post] = isolator_make_input(P)
%ISOLATOR_MAKE_INPUT  Capsule accelerometer input shared by both models.
% Builds the synthetic pulse (or loads measured data) once, from the inputs
% in isolator_inputs.m, so the two model scripts are driven by the
% identical sampled record. Output a is the capsule accelerometer reading
% (specific force, +1 g at rest) in m/s^2 on a uniform grid t.
%
% The grid resolves the fastest system in the study: the stiffest swept
% isolator on the lightest mass AND the Brinkley body model (wn ~ 60 rad/s,
% about 10 Hz, which is faster than any isolator in the sweep).
g = P.g;
fn_iso = sqrt(max(P.k_vec)/min(P.m_list))/(2*pi);
fn_br  = 0;
if P.use_brinkley, fn_br = P.BTAB.(P.brinkley_axis){1}/(2*pi); end
fn_max = max(fn_iso, fn_br);

if strcmp(P.pulse, 'data')
    [t, a, t_imp, dt] = load_accel(P.data_file, P.data_units, P.data_sign, P.t_imp_data, ...
                                   P.detect_g, P.pad_post, P.t_post, fn_max, P.ppc, g);
    T = NaN;
else
    [t, a, T, t_imp, dt] = make_pulse(P.pulse, P.Vi, P.a_in_pk, P.t_pre, P.t_post, ...
                                      fn_max, P.n_min, P.ppc, g);
end
a0   = mean(a(t < t_imp));      % pre-impact reading (equilibrium)
post = t >= t_imp;
end

function [t, a, t_imp, dt] = load_accel(file, units, sgn, t_imp, detect_g, ...
                                         pad_post, t_post, fn, ppc, g)
% Loads measured capsule accelerometer data, cleans it, finds impact, and
% resamples it onto a uniform grid fine enough for the stiffest swept system.
    d = readmatrix(file);
    d = d(all(isfinite(d(:,1:2)), 2), 1:2);
    d = sortrows(d, 1);
    [tr, iu] = unique(d(:,1));
    ar = sgn * d(iu, 2);
    switch units
        case 'g',     ar = ar * g;
        case 'm/s^2'
        otherwise, error('data_units must be ''g'' or ''m/s^2''');
    end
    tr = tr - tr(1);

    % Sanity check on sign and units: pre-impact should read about +1 g
    nb   = max(20, round(0.05*numel(ar)));
    base = median(ar(1:nb));
    if abs(base/g - 1) > 0.25
        warning(['Start of record reads %.2f g, not about +1 g. ' ...
                 'Check data_units, data_sign, and that the file starts before impact.'], base/g);
    end

    % Impact time
    if isempty(t_imp)
        dev = abs(ar - base);
        i1  = find(dev > detect_g*g, 1);
        if isempty(i1), error('No impact found: nothing deviates %g g from baseline.', detect_g); end
        while i1 > 1 && dev(i1-1) > 0.1*detect_g*g, i1 = i1 - 1; end   % back up to onset
        t_imp = tr(i1);
    else
        t_imp = t_imp - d(1,1);          % same shift as the time vector
    end
    if t_imp < 0.02
        warning('Less than 20 ms of data before impact; the equilibrium estimate may be poor.');
    end

    % Optional padding with the floating (post-impact) reading
    if pad_post && tr(end) < t_imp + t_post
        nt    = max(20, round(0.05*numel(ar)));
        a_end = mean(ar(end-nt+1:end));
        dt_r  = median(diff(tr));
        tpad  = (tr(end)+dt_r : dt_r : t_imp + t_post)';
        tr    = [tr; tpad];
        ar    = [ar; a_end*ones(size(tpad))];
        fprintf('Padded record with %.2f g out to %.2f s\n', a_end/g, tr(end));
    end

    % Uniform resampling
    dt_data = median(diff(tr));
    dt      = min(dt_data, 1/(ppc*fn));
    if dt < dt_data
        fprintf(['Data sampled at %.0f Hz; interpolating to %.0f Hz to resolve the ' ...
                 'stiffest swept system (no new content is added).\n'], 1/dt_data, 1/dt);
    end
    t = (0:dt:tr(end))';
    a = interp1(tr, ar, t, 'linear');
end

function [t, a, T, t_imp, dt] = make_pulse(type, Vi, a_pk, t_pre, t_post, fn, n_min, ppc, g)
% Capsule accelerometer input: 1 g, pulse, then 1 g.
% Pulse area above 1 g equals Vi.
    A = (a_pk - 1)*g;
    if A <= 0, error('a_in_pk must be greater than 1 g'); end
    switch type
        case 'halfsine', T = pi*Vi/(2*A);
        case 'rect',     T = Vi/A;
        otherwise, error('pulse must be ''halfsine'' or ''rect''');
    end

    N  = ceil(max(n_min, ppc*fn*T));  % resolve pulse and fastest system
    dt = T/N;

    n_pre  = round(t_pre/dt);
    n_post = round(t_post/dt);
    t      = (0:(n_pre + N + n_post))' * dt;
    t_imp  = n_pre*dt;

    a = g*ones(size(t));
    switch type
        case 'halfsine'
            idx = n_pre + (0:N) + 1;
            a(idx) = g + A*sin(pi*(0:N)'*dt/T);
        case 'rect'
            a(n_pre + (1:N)) = g + A;
    end
end