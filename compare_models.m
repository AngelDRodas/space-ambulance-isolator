%% compare_models.m
% Cross-validation of the two isolator models.
%   1. Run ShockAttentuationSpaceAmbulanceBrinkley.m   -> results_vibrations.mat
%   2. Run ShockAttentuationSpaceAmbulanceBondGraphs.m -> results_bondgraph.mat
%   3. Run this script.
%
% Both models solve the same physics (base-excited spring-damper + Brinkley
% body model) by different routes: an exact discrete state-space filter
% (vibrations) and RK4 integration of bond-graph state equations. If both
% are right, every metric agrees to within the time-step error, and that
% error falls about 4x when the time step is halved (both are 2nd-order
% accurate in how they treat the sampled input).
%
% Reports, per mass:
%   - relative error of stroke, peak accel, impulse and beta over the whole
%     (k, c) sweep: max |e|, RMS e, and where the max occurs
%   - pass/fail agreement (designs classed differently should sit on a limit)
%   - whether both models select the same design
%   - time-history error at the check design (k_check, c_check):
%     normalised max error  max|x_bg - x_vib| / max|x_vib|  and NRMSE
% Writes model_comparison.csv and model_comparison_table.tex.

clear; clc; close all
V = load('results_vibrations.mat');
B = load('results_bondgraph.mat');

% ---- same case? ----
assert(isequal(V.m_list, B.m_list) && isequal(V.k_vec, B.k_vec) && ...
       isequal(V.c_vec, B.c_vec), 'The two runs used different sweeps.');
assert(isequal(size(V.t), size(B.t)) && max(abs(V.a - B.a)) == 0, ...
       'The two runs used different input records. Rerun both scripts.');
assert(V.chk.k == B.chk.k && V.chk.c == B.chk.c, ...
       ['The two runs saved different check designs (%.0f/%.0f vs %.0f/%.0f). ' ...
        'Set P.k_check and P.c_check in isolator_inputs.m.'], V.chk.k, V.chk.c, B.chk.k, B.chk.c);
fprintf('Input: %s, Vi = %.2f m/s, dt = %.3f ms, %d designs x %d masses\n\n', ...
        V.IN.pulse, V.IN.Vi, 1e3*V.dt, numel(V.k_vec)*numel(V.c_vec), numel(V.m_list));

names = {'stroke', 'accel', 'impulse', 'beta'};
units = {'m', 'g', 'g*s', '-'};
VV = {V.STROKE, V.ACC, V.IMP, V.BETA};
BB = {B.STROKE, B.ACC, B.IMP, B.BETA};
[KK, CC] = meshgrid(V.k_vec, V.c_vec);
nm = numel(V.m_list);

rows = {};   % for CSV / LaTeX: mass, quantity, max|e|, rms e, at k, at c
fprintf('Sweep-wide relative error, (bond graph - vibrations) / vibrations\n');
fprintf('%8s %9s %12s %12s %10s %10s\n', 'm [kg]', 'metric', 'max|e| [%]', ...
        'RMS e [%]', 'at k', 'at c');
for im = 1:nm
    for q = 1:numel(names)
        x = VV{q}(:,:,im);  y = BB{q}(:,:,im);
        if all(x(:) == 0), continue; end          % metric switched off
        e = (y - x) ./ x;
        [emax, i] = max(abs(e(:)));
        erms = sqrt(mean(e(:).^2));
        fprintf('%8.2f %9s %12.2e %12.2e %10.0f %10.0f\n', V.m_list(im), names{q}, ...
                100*emax, 100*erms, KK(i), CC(i));
        rows(end+1, :) = {V.m_list(im), names{q}, 100*emax, 100*erms, KK(i), CC(i)}; %#ok<SAGROW>
    end
end

% ---- pass/fail agreement ----
fprintf('\nPass/fail agreement\n');
for im = 1:nm
    d = xor(V.PASS(:,:,im), B.PASS(:,:,im));
    fprintf('  m = %.2f kg: %d of %d designs classed differently\n', ...
            V.m_list(im), nnz(d), numel(d));
    if any(d(:))
        % how close were those designs to a limit? (worst ratio, vibrations)
        ratio = worst_ratio(V, im);
        fprintf('     their worst metric/limit ratio spans %.5f - %.5f (1 = on the limit)\n', ...
                min(ratio(d)), max(ratio(d)));
    end
end
dall = xor(all(V.PASS, 3), all(B.PASS, 3));
fprintf('  all masses: %d designs classed differently\n', nnz(dall));
fprintf('\nSelected design: vibrations k = %.0f, c = %.0f | bond graph k = %.0f, c = %.0f  -> %s\n', ...
        V.k_best, V.c_best, B.k_best, B.c_best, ...
        tern(V.k_best == B.k_best && V.c_best == B.c_best, 'SAME', 'DIFFERENT (check near-ties)'));

% ---- check-design time histories ----
fprintf('\nTime histories at check design k = %.0f N/m, c = %.0f N*s/m\n', V.chk.k, V.chk.c);
fprintf('%8s %9s %14s %12s\n', 'm [kg]', 'signal', 'norm. max e', 'NRMSE');
sig  = {'stroke', 'am', 'DR'};
slab = {'stroke', 'accel', 'DR'};
for im = 1:nm
    for s = 1:numel(sig)
        x = V.chk.(sig{s})(:,im);  y = B.chk.(sig{s})(:,im);
        nmax = max(abs(y - x)) / max(abs(x));
        nrms = sqrt(mean((y - x).^2)) / (max(x) - min(x));
        fprintf('%8.2f %9s %14.2e %12.2e\n', V.m_list(im), slab{s}, nmax, nrms);
        rows(end+1, :) = {V.m_list(im), ['history ' slab{s}], 100*nmax, 100*nrms, V.chk.k, V.chk.c}; %#ok<SAGROW>
    end
end

% check-design performance (from the vibrations histories; bond graph agrees)
lim = V.IN.BTAB.(V.IN.brinkley_axis){3 + strcmp(V.IN.crew_condition, 'ill')};
fprintf('\nCheck design performance (vibrations model)\n');
fprintf('%8s %10s %10s %8s\n', 'm [kg]', 'stroke[m]', 'accel[g]', 'beta');
for im = 1:nm
    DR = V.chk.DR(:,im);
    fprintf('%8.2f %10.3f %10.2f %8.3f\n', V.m_list(im), max(V.chk.stroke(:,im)), ...
            max(abs(V.chk.am(:,im)))/V.IN.g, max([max(DR)/lim(1), min(DR)/lim(2), 0]));
end

% ---- CSV ----
fid = fopen('model_comparison.csv', 'w');
fprintf(fid, 'mass_kg,quantity,max_abs_err_pct,rms_err_pct,at_k_Npm,at_c_Nspm\n');
for r = 1:size(rows, 1)
    fprintf(fid, '%.2f,%s,%.4e,%.4e,%.0f,%.0f\n', rows{r, :});
end
fclose(fid);

% ---- LaTeX table (sweep-wide metric errors) ----
fid = fopen('model_comparison_table.tex', 'w');
fprintf(fid, '%% Generated by compare_models.m\n');
fprintf(fid, '\\begin{table}[htbp]\n  \\centering\n');
fprintf(fid, ['  \\caption{Cross-validation of the vibrations (exact discrete filter) and ' ...
              'bond-graph (RK4) isolator models over the full $k$--$c$ sweep ' ...
              '(%d designs per mass, %s pulse, $V_i = %.1f$~m/s, ' ...
              '$\\Delta t = %.2f$~ms). Error is (bond graph $-$ vibrations)/vibrations.}\n'], ...
        numel(V.k_vec)*numel(V.c_vec), strrep(V.IN.pulse, 'halfsine', 'half-sine'), ...
        V.IN.Vi, 1e3*V.dt);
fprintf(fid, '  \\label{tab:transport_model_cross_validation}\n');
fprintf(fid, '  \\begin{tabular}{llrr}\n    \\hline\n');
fprintf(fid, '    Isolated mass (kg) & Metric & Max $|e|$ (\\%%) & RMS $e$ (\\%%) \\\\\n    \\hline\n');
for r = 1:size(rows, 1)
    if strncmp(rows{r, 2}, 'history', 7), continue; end
    fprintf(fid, '    %.1f & %s & %.1e & %.1e \\\\\n', rows{r, 1}, rows{r, 2}, rows{r, 3}, rows{r, 4});
end
fprintf(fid, '    \\hline\n  \\end{tabular}\n\\end{table}\n');
fclose(fid);
fprintf('\nWrote model_comparison.csv and model_comparison_table.tex\n');

% ---- plots ----
if usejava('desktop') || ~isempty(getenv('DISPLAY')) || ispc || ismac
    for im = 1:nm
        figure('Name', sprintf('Model error, m = %.1f kg', V.m_list(im)), 'Color', 'w');
        for q = 1:numel(names)
            x = VV{q}(:,:,im);  y = BB{q}(:,:,im);
            subplot(2,2,q);
            if all(x(:) == 0), axis off; title([names{q} ' (off)']); continue; end
            imagesc(V.k_vec, V.c_vec, log10(abs((y - x)./x) + eps)); axis xy; colorbar
            xlabel('k [N/m]'); ylabel('c [N\cdot{}s/m]');
            title(sprintf('log_{10}|rel. error|, %s', names{q}));
        end
    end
    figure('Name', 'Check design time histories', 'Color', 'w');
    for s = 1:numel(sig)
        subplot(numel(sig), 2, 2*s-1); hold on
        plot(V.t, V.chk.(sig{s}), '-', 'LineWidth', 1.4);
        set(gca, 'ColorOrderIndex', 1);
        plot(B.t, B.chk.(sig{s}), '--', 'LineWidth', 1.4);
        grid on; ylabel(slab{s});
        if s == 1, title('solid: vibrations, dashed: bond graph'); end
        subplot(numel(sig), 2, 2*s); plot(V.t, B.chk.(sig{s}) - V.chk.(sig{s}));
        grid on; ylabel(['\Delta ' slab{s}]);
        if s == 1, title('bond graph - vibrations'); end
    end
    xlabel('t [s]');
end

function r = worst_ratio(S, im)
% worst metric/limit ratio for one mass (hard limits, legacy limits if on)
    IN = S.IN;
    r = S.STROKE(:,:,im) / IN.stroke_max;
    if IN.use_legacy
        r = max(r, S.ACC(:,:,im) / IN.accel_max);
        r = max(r, S.IMP(:,:,im) / IN.impulse_max);
    end
    if IN.use_brinkley, r = max(r, S.BETA(:,:,im)); end
end

function o = tern(c, a, b)
    if c, o = a; else, o = b; end
end
