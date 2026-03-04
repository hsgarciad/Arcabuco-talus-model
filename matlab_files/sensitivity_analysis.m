clear; clc;

C10_measured   = 1.89e6;
C26_measured   = 1.14e7;
ratio_measured = C26_measured / C10_measured;

%% Field constraints — update these with your site data
field.H_min       = 2;     % m
field.H_max       = 50;    % m — from geological mapping
field.L_cliff_min = 10;    % m — update from DEM/field
field.L_cliff_max = 100;   % m
field.theta_min   = 8;     % degrees
field.theta_max   = 25;    % degrees — from longitudinal profiles
field.S_min       = 0.5;   % ha
field.S_max       = 100;   % ha
field.rp_min      = 10;    % yr
field.rp_max      = 1000;  % yr

%% Search bounds
% [log10(e_back), log10(S), L_cliff, theta, rp]
bounds_min = [log10(1e-7), log10(field.S_min), field.L_cliff_min, ...
              field.theta_min, field.rp_min];
bounds_max = [log10(1e-3), log10(field.S_max), field.L_cliff_max, ...
              field.theta_max, field.rp_max];

%% Misfit function using reduced tsim
    function misfit = fast_misfit(p, C10_meas, C26_meas, ratio_meas, field)
        e_back  = 10^p(1);
        S       = 10^p(2);
        L_cliff = p(3);
        theta   = p(4);
        rp      = max(1, round(p(5)));

        % Hard field bounds
        if S       < field.S_min       || S       > field.S_max       || ...
           L_cliff < field.L_cliff_min || L_cliff > field.L_cliff_max || ...
           theta   < field.theta_min   || theta   > field.theta_max   || ...
           rp      < field.rp_min      || rp      > field.rp_max      || ...
           e_back  <= 0
            misfit = 1e6;
            return
        end

        base.e_back  = e_back;
        base.S       = S;
        base.L_cliff = L_cliff;
        base.theta   = theta;
        base.rp      = rp;
        base.Ast     = 0.03;
        base.tsim    = 5e5;   % 500 kyr — fast during search

        try
            [C10, C26, H] = RockFall_v3_sim(base);
        catch
            misfit = 1e6;
            return
        end

        % Hard pile height constraint
        if H < field.H_min || H > field.H_max
            misfit = 1e6;
            return
        end

        ratio_model = C26 / C10;
        res10   = (C10         - C10_meas)   / C10_meas;
        res26   = (C26         - C26_meas)   / C26_meas;
        res_rat = (ratio_model - ratio_meas) / ratio_meas;
        misfit  = res10^2 + res26^2 + res_rat^2;
    end

%% Multi-start optimization
n_starts     = 100;
rng(42);
all_solutions = [];

options = optimset('TolX', 1e-4, 'TolFun', 1e-6, 'Display', 'off', ...
                   'MaxIter', 500, 'MaxFunEvals', 2000);

fprintf('Running %d optimizations (fast mode: tsim=500kyr)...\n\n', n_starts);
fprintf('%-6s %-12s %-8s %-10s %-8s %-6s %-10s %-10s %-10s %-10s %-8s\n', ...
    'Run', 'e_back', 'S(ha)', 'L_cliff', 'theta', 'rp', ...
    'C10 err%', 'C26 err%', 'Ratio err%', 'Misfit', 'H(m)');

for k = 1:n_starts

    % First start from sensitivity best values
    if k == 1
        p0 = [log10(4.833e-6), log10(1.526), 26.32, 8.0, 89];
    else
        p0 = bounds_min + rand(1,5) .* (bounds_max - bounds_min);
    end

    try
        obj   = @(p) fast_misfit(p, C10_measured, C26_measured, ...
                                 ratio_measured, field);
        p_opt = fminsearch(obj, p0, options);

        % Extract
        e_back_sol  = 10^p_opt(1);
        S_sol       = 10^p_opt(2);
        L_cliff_sol = p_opt(3);
        theta_sol   = p_opt(4);
        rp_sol      = max(1, round(p_opt(5)));

        % Evaluate with full tsim for accurate final values
        base.e_back  = e_back_sol;
        base.S       = S_sol;
        base.L_cliff = L_cliff_sol;
        base.theta   = theta_sol;
        base.rp      = rp_sol;
        base.Ast     = 0.03;
        base.tsim    = 1e6;   % full run for reporting

        [C10, C26, H] = RockFall_v3_eback(base);
        ratio_model   = C26 / C10;

        res10   = (C10         - C10_measured)   / C10_measured;
        res26   = (C26         - C26_measured)   / C26_measured;
        res_rat = (ratio_model - ratio_measured) / ratio_measured;
        misfit  = res10^2 + res26^2 + res_rat^2;

        err10   = abs(res10)   * 100;
        err26   = abs(res26)   * 100;
        err_rat = abs(res_rat) * 100;

        fprintf('%-6d %-12.3e %-8.3f %-10.1f %-8.1f %-6d %-10.1f %-10.1f %-10.1f %-10.4f %-8.1f\n', ...
            k, e_back_sol, S_sol, L_cliff_sol, theta_sol, rp_sol, ...
            err10, err26, err_rat, misfit, H);

        all_solutions(end+1,:) = [e_back_sol, S_sol, L_cliff_sol, ...
                                  theta_sol, rp_sol, misfit, ...
                                  err10, err26, err_rat, H];
    catch
        continue
    end
end

if isempty(all_solutions)
    fprintf('\nNo solutions converged. Check run_rockfall.m accepts tsim from params.\n');
    return
end

%% Sort by misfit
[~, sort_idx] = sort(all_solutions(:,6));
all_solutions = all_solutions(sort_idx,:);

%% Statistics
fprintf('\n=== MISFIT STATISTICS ===\n');
fprintf('  Minimum : %.4f\n', min(all_solutions(:,6)));
fprintf('  Median  : %.4f\n', median(all_solutions(:,6)));
fprintf('  Maximum : %.4f\n', max(all_solutions(:,6)));

%% Top 10
fprintf('\n=== TOP 10 SOLUTIONS ===\n');
fprintf('%-6s %-12s %-8s %-10s %-8s %-6s %-10s %-10s %-10s %-10s %-8s\n', ...
    'Rank', 'e_back', 'S(ha)', 'L_cliff', 'theta', 'rp', ...
    'C10 err%', 'C26 err%', 'Ratio err%', 'Misfit', 'H(m)');

for k = 1:min(10, size(all_solutions,1))
    fprintf('%-6d %-12.3e %-8.3f %-10.1f %-8.1f %-6d %-10.1f %-10.1f %-10.1f %-10.4f %-8.1f\n', ...
        k, all_solutions(k,1), all_solutions(k,2), all_solutions(k,3), ...
        all_solutions(k,4), round(all_solutions(k,5)), ...
        all_solutions(k,7), all_solutions(k,8), all_solutions(k,9), ...
        all_solutions(k,6), all_solutions(k,10));
end

%% Accept solutions within 3x best misfit
best_misfit = min(all_solutions(:,6));
threshold   = best_misfit * 3;
accepted    = all_solutions(all_solutions(:,6) <= threshold, :);

fprintf('\n=== ACCEPTED SOLUTIONS (within 3x best misfit = %.4f) ===\n', threshold);
fprintf('Found %d solutions\n\n', size(accepted,1));

%% Parameter ranges across accepted solutions
fprintf('%-14s %12s %12s %12s %12s\n', 'Parameter', 'Min', 'Max', 'Mean', 'Std');
param_labels = {'e_back', 'S (ha)', 'L_cliff (m)', 'theta (deg)', 'rp (yr)'};
for p = 1:5
    vals = accepted(:,p);
    fprintf('%-14s %12.3e %12.3e %12.3e %12.3e\n', ...
        param_labels{p}, min(vals), max(vals), mean(vals), std(vals));
end

%% Scatter plot matrix
if size(accepted,1) > 1
    figure('Position', [50 50 1200 1000]);
    sol_plot      = accepted(:,1:5);
    sol_plot(:,1) = log10(sol_plot(:,1));
    sol_plot(:,2) = log10(sol_plot(:,2));
    param_short   = {'log(e\_back)', 'log(S)', 'L\_cliff', 'theta', 'rp'};

    plot_idx = 1;
    for i = 1:5
        for j = i+1:5
            subplot(4, 5, plot_idx);
            scatter(sol_plot(:,j), sol_plot(:,i), 60, accepted(:,6), 'filled');
            colormap(flipud(hot)); colorbar;
            xlabel(param_short{j}, 'Interpreter', 'tex');
            ylabel(param_short{i}, 'Interpreter', 'tex');
            grid on;
            plot_idx = plot_idx + 1;
        end
    end
    sgtitle('Solution cloud — color = misfit (darker = better)');
end

%% Misfit histogram
figure('Position', [50 50 500 350]);
histogram(all_solutions(:,6), 20, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor','w');
xline(best_misfit, 'r-', 'Best', 'LineWidth', 2);
xline(threshold,   'g--', '3x Best', 'LineWidth', 1.5);
xlabel('Misfit'); ylabel('Count');
title('Misfit distribution across all starts');
grid on;