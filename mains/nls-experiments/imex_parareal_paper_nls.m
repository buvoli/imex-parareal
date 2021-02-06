%% -- Include stepper files --------------------------------------------------------------------------------------------
addpath(genpath(fullfile(pwd(), '../../stepper/common')));
addpath(fullfile(pwd(), '../../stepper/rk'));
addpath(fullfile(pwd(), '../../stepper/parareal'));

%% -- Load Equation ----------------------------------------------------------------------------------------------------
equation_name = 'nls';

    cwd = pwd();
    cd(['../../stepper/equations/',equation_name]);
    init
    cd(cwd)

%% -- Reference Solution -----------------------------------------------------------------------------------------------
erk_options = struct('coeffGenerator', @IMRK4, 'parameters', pars,'max_ts_to_store', 2);
Nt_ref = 2^19;
[TS, YS, tcpu_erkG, tccpu_erkG] = IMRK(LF, NF, tspan, y0, Nt_ref, erk_options); 
y_ref = YS(:, end);
    
%% -- Experiment 1: Vary Nt, Fix Nt ------------------------------------------------------------------------------------ 
steps    = 2 .^ (6:18);
ks       = 0:6;
Nts      = [512 1024 2048 4096 8192];
num_fine = 16;

% init structs
results_parareal(length(Nts), length(steps), length(ks)).speedup = 0;   [results_parareal.speedup] = deal(0);       
results_parareal(length(Nts), length(steps), length(ks)).error = NaN;   [results_parareal.error] = deal(NaN);
results_parareal(length(Nts), length(steps), length(ks)).optimal_time = NaN; [results_parareal.optimal_time] = deal(NaN);

results_rk(length(steps)).error = NaN; [results_rk.error] = deal(NaN);
result_rk(length(steps)).time = 0; [results_rk.time] = deal(0);

%% -- run fine integrator ----------------------------------------------------------------------------------------------
for j = 1 : length(steps)
    fprintf('IMEX-RK4 (Ns = %i)\n', steps(j));
    [t_rkf, y_rkf, tcpu] = IMRK(LF, NF,tspan, y0, steps(j), struct('parameters', pars, 'coeffGenerator', @IMRK4, 'max_ts_to_store', 2)); 
    error = error_filters{1}(y_rkf(:, end), y_ref);
     
    results_rk(j).error = error;
    results_rk(j).time  = tcpu;
end

%% -- run parareal -----------------------------------------------------------------------------------------------------
for i = 1 : length(Nts)
    
    [~, first_step_index] = find(steps >= Nts(i));
    
    for j = first_step_index : length(steps)
        
        for k = 1 : length(ks)
        
            fprintf('IMEX-PR (Nt = %i, k = %i, Ns = %i)\n', Nts(i), ks(k), steps(j));

            pr_options = struct(                        ...
                'coarse',       IMEXRKWorker(@IMRK3),   ...
                'fine',         IMEXRKWorker(@IMRK4),   ...
                'num_iter',     ks(k),                  ... 
                'num_coarse',   1,                      ...
                'num_fine',     num_fine,               ...
                'num_proc',     Nts(i) / num_fine       ...
            );

            Nb = steps(j) / (pr_options.num_proc * pr_options.num_fine);
            
            % -- run parareal ------------------------------------------------------------------------------------------
            PR = parareal(pr_options);
            [t_pr, y_pr] = PR.solve(LF, @(t,y) NF(t, y, pars), tspan, y0, Nb);        
            results_parareal(i,j,k).error = error_filters{1}(y_pr, y_ref);
            results_parareal(i,j,k).speedup = pararealSpeedup(pr_options);
            results_parareal(i,j,k).optimal_time = results_rk(j).time / results_parareal(i,j,k).speedup;
            
        end
        
    end

end

%% -- save data --------------------------------------------------------------------------------------------------------
basename = ['imex-nls-paper-experiment-A-', num2str(floor(posixtime(datetime('now'))))];
save(fullfile('data', [basename, '.mat']))

%% -- plot settings -------------------------------------------------------
marker_styles = {'o', '+', '*', 'x', 's', 'd', '^'};
markerStyle   = @(j) marker_styles{mod(j-1, length(marker_styles)) + 1};
marker_size   = 10;

line_styles   = {'-',':','-.','--','--*','--+', '--d'};
lineStyle   = @(j) line_styles{mod(j-1, length(line_styles)) + 1};

%% -- generate plots ---------------------------------------------------------------------------------------------------
hs = diff(tspan) ./ steps;
scale_times = true;

if(scale_times) % scale times relative to first non NaN fine rk
    [~, ind] = find(~isnan([results_rk.error]), 1);
    tsf = results_rk(ind).time; 
else
    tsf = 1;
end

legend_cell = [ {'IMEX-RK4'}, arrayfun(@(k) ['IMEX-PR (k=', num2str(k), ')'], ks, 'UniformOutput', false) ];
rk_errors = [results_rk.error];
rk_times = [results_rk.time];

for i = 1 : length(Nts)
    fh = figure(1);
    
    title_str = sprintf('{\\boldmath $N_T=%i$, $N_F=%i$, $N_P=%i$}', Nts(i), num_fine, Nts(i) / num_fine);
    p_error_matrix   = reshape([results_parareal(i,:,:).error], [length(steps), length(ks)]);
    p_times_matrix   = reshape([results_parareal(i,:,:).optimal_time], [length(steps), length(ks)]);
    
    % -- accuracy plots ------------------------------------------------------------------------------------------------
    basename = ['imex-nls-accuracy-NT', num2str(Nts(i))];
    loglog(hs, rk_errors, '.-k', 'linewidth', 3, 'MarkerSize', 30); hold on;
    for j = 1 : size(p_error_matrix,2)
        p1 = loglog(hs, p_error_matrix(:,j), ['-', markerStyle(j)],  'linewidth', 3, 'MarkerSize', marker_size);
        set(p1, 'MarkerFaceColor', get(p1, 'color'));
    end
    hold off;
    if(i == 1)
        legend(legend_cell, 'Location', 'southeast'); legend box off;
    end
    axis([min(hs) max(hs) 1e-9 1e1]);
    xlabel('stepsize $\Delta t$', 'interpreter', 'latex'); ylabel('relative error', 'interpreter', 'latex');
    title(title_str, 'interpreter', 'Latex');
    set(gca, 'FontSize', 15, 'FontName', 'Minion Pro');
    exportFigure(fh, struct('SavePath', fullfile('figures', basename), 'Format', 'jpg', 'PaperPosition', [0 0 12, 14]))
    
    % -- efficiency plots ----------------------------------------------------------------------------------------------
    basename = ['imex-nls-efficiency-NT', num2str(Nts(i))];
    loglog(rk_times / tsf, rk_errors, '.-k', 'linewidth', 3, 'MarkerSize', 30); hold on;
    for j = 1 : size(p_error_matrix,2)
        loglog(p_times_matrix(:,j) / tsf, p_error_matrix(:,j), lineStyle(j), 'linewidth', 3);
    end
    hold off;
    legend(legend_cell, 'Location', 'northeast'); %legend box off;
    if(scale_times)
        xlabel('relative time');
        axis([1 max(rk_times / tsf) 1e-9 1e1]);
    else
        xlabel('time (sec)');
        axis([min(rk_times) max(rk_times)/sf 1e-9 1e1]);
    end
    ylabel('relative error');    
    title(title_str, 'interpreter', 'Latex');
    set(gca, 'FontSize', 17, 'FontName', 'Minion Pro');
    exportFigure(fh, struct('SavePath', fullfile('figures', basename), 'Format', 'jpg', 'PaperPosition', [0 0 16, 14]))

end   

%% -- generate adaptive efficiency diagrams ----------------------------------------------------------------------------
libpfasst_data = csvread('libpfasst-data-no-headers.csv');
normalization  = libpfasst_data(2,2); % normalize all times relative to first convergent timestep of fine integrator

rk_errors = [results_rk.error];
rk_times  = libpfasst_data(:,2) / normalization;

% -- Nt = 2048, k=3 ----------------------------------------------------------------------------------------------------
p_k3_error = [results_parareal(3, :, 4).error];
p_k3_times = libpfasst_data(:,3) / normalization;
p_k3_times_theoretical = rk_times / 16.176935229067933;

% -- Nt = 2048, adaptive k<=3 ------------------------------------------------------------------------------------------
p_ka_error = libpfasst_data(:, 5);
p_ka_times = libpfasst_data(:, 4) / normalization;
p_ka_times_theoretical = rk_times(:) ./ libpfasst_data(:,6);

fh = figure(1); clf;
loglog(rk_times, rk_errors, '.-k', 'linewidth', 3, 'MarkerSize', 30); hold on;
loglog(p_ka_times, p_ka_error, 's-', 'color', [0.635000000000000   0.078000000000000   0.184000000000000], 'MarkerFaceColor', [0.635000000000000   0.078000000000000   0.184000000000000], 'linewidth', 3, 'MarkerSize', 13);
loglog(p_k3_times, p_k3_error, 'd-', 'color', [0.466000000000000   0.674000000000000   0.188000000000000], 'MarkerFaceColor', [0.466000000000000   0.674000000000000   0.188000000000000], 'linewidth', 3, 'MarkerSize', 7);

loglog(p_ka_times_theoretical, p_ka_error, 's-.', 'color', .3 * [1 1 1], 'MarkerFaceColor', .3 * [1 1 1], 'linewidth', 3,  'MarkerSize', 13);
loglog(p_k3_times_theoretical, p_k3_error, 'd--', 'color', .6 * [1 1 1], 'MarkerFaceColor', .6 * [1 1 1], 'linewidth', 3,  'MarkerSize', 7);

title_str = sprintf('{\\boldmath $N_T=%i$, $N_F=%i$, $N_P=%i$}', Nts(3), num_fine, Nts(3) / num_fine);
title(title_str, 'interpreter', 'Latex');
legend({'IMEX-RK4', 'IMEX-PR, $k\le3$', 'IMEX-PR, $k=3$', 'IMEX-PR, $k\le3$ (theory)', 'IMEX-PR, $k=3$ (theory)'}, 'interpreter', 'latex', 'Location', 'NorthEast'); legend box off;
set(gca, 'FontSize', 17, 'FontName', 'Minion Pro');
xlabel('relative time'); 
ylabel('relative error');    
axis([1 rk_times(end) 1e-9 1e1])

exportFigure(fh, struct('SavePath', fullfile('figures', 'imex-nls-efficiency-adaptive'), 'Format', 'jpg', 'PaperPosition', [0 0 16, 14]))

%% -- Experiment 2: Fix Nt, vary Nf ------------------------------------------------------------------------------------ 
Nt_fixed  = 512;
k_fixed   = 3;
Nfs = [4 8 16 32];

% init structs
resultsB_parareal(length(steps), length(Nfs)).speedup = 0;   [resultsB_parareal.speedup] = deal(0);       
resultsB_parareal(length(steps), length(Nfs)).error = NaN;   [resultsB_parareal.error] = deal(NaN);
resultsB_parareal(length(steps), length(Nfs)).optimal_time = NaN; [resultsB_parareal.optimal_time] = deal(NaN);

%% -- run parareal -----------------------------------------------------------------------------------------------------
[~, first_step_index] = find(steps >= Nt_fixed);
for j = first_step_index : length(steps)
    
    for k = 1 : length(Nfs)
        
        fprintf('IMEX-PR (Nt = %i, Nf = %i, k = %i, Np = %i)\n', Nt_fixed, k_fixed, steps(j), Nt_fixed / Nfs(k));

        pr_options = struct(                        ...
            'coarse',       IMEXRKWorker(@IMRK3),   ...
            'fine',         IMEXRKWorker(@IMRK4),   ...
            'num_iter',     k_fixed,                ... 
            'num_coarse',   1,                      ...
            'num_fine',     Nfs(k),                 ...
            'num_proc',     Nt_fixed / Nfs(k)       ...
        );

        Nb = steps(j) / (pr_options.num_proc * pr_options.num_fine);

        % -- run parareal ----------------------------------------------------------------------------------------------
        PR = parareal(pr_options);
        [t_pr, y_pr] = PR.solve(LF, @(t,y) NF(t, y, pars), tspan, y0, Nb);        
        resultsB_parareal(j,k).error = error_filters{1}(y_pr, y_ref);
        resultsB_parareal(j,k).speedup = pararealSpeedup(pr_options);
        resultsB_parareal(j,k).optimal_time = results_rk(j).time / resultsB_parareal(j,k).speedup;
    
    end
            
end

%% -- save data --------------------------------------------------------------------------------------------------------
basename = ['imex-nls-paper-experiment-B-', num2str(floor(posixtime(datetime('now'))))];
save(fullfile('data', [basename, '.mat']), 'Nt_fixed', 'k_fixed', 'Nfs', 'resultsB_parareal', 'steps');

%% -- generate plots ---------------------------------------------------------------------------------------------------
hs = diff(tspan) ./ steps;
%legend_cell = [ {'IMEX-RK4'}, arrayfun(@(k) ['IMEX-PR (N_F=', num2str(k), ')'], Nfs, 'UniformOutput', false) ];
legend_cell = [ {'IMEX-RK4'}, arrayfun(@(k) ['IMEX-PR (N_P=', num2str(k), ')'], Nt_fixed ./ Nfs, 'UniformOutput', false) ];
rk_errors = [results_rk.error];
rk_times = [results_rk.time];

fh = figure(1);

title_str = sprintf('{\\boldmath $N_T=%i$, $k=%i$}', Nt_fixed, k_fixed);
p_error_matrix   = reshape([resultsB_parareal.error], [length(steps), length(Nfs)]);
p_times_matrix   = reshape([resultsB_parareal.optimal_time], [length(steps), length(Nfs)]);

% -- accuracy plots ----------------------------------------------------------------------------------------------------
basename = ['imex-nls-NFS-accuracy-NT', num2str(Nt_fixed)];
loglog(hs, rk_errors, '.-k', 'linewidth', 3, 'MarkerSize', 30); hold on;
for i = 1 : size(p_error_matrix, 2)
    loglog(hs, p_error_matrix(:,i), [markerStyle(i), '-'], 'linewidth', 3, 'MarkerSize', marker_size);
end
hold off;
legend(legend_cell, 'Location', 'southeast'); legend box off;
axis([min(hs) max(hs) 1e-9 1e1]);
xlabel('stepsize $\Delta t$', 'interpreter', 'latex'); ylabel('relative error', 'interpreter', 'latex');
title(title_str, 'interpreter', 'Latex');
set(gca, 'FontSize', 15, 'FontName', 'Minion Pro');
exportFigure(fh, struct('SavePath', fullfile('figures', basename), 'Format', 'jpg', 'PaperPosition', [0 0 12, 14]))

% -- efficiency plots --------------------------------------------------------------------------------------------------
clf;
basename = ['imex-nls-NFS-efficiency-NT', num2str(Nt_fixed)];
loglog(rk_times, rk_errors, '.-k', 'linewidth', 3, 'MarkerSize', 30); hold on;
for i = 1 : size(p_error_matrix, 2)
    p1 = loglog(p_times_matrix(:,i), p_error_matrix(:,i), [markerStyle(i), '--'], 'linewidth', 3, 'MarkerSize', marker_size);
end
hold off;
legend(legend_cell, 'Location', 'southwest'); legend box off;
axis([min(rk_times) max(rk_times) 1e-9 1e1]);
title(title_str, 'interpreter', 'Latex');
xlabel('time (sec)'); ylabel('relative error');
set(gca, 'FontSize', 15, 'FontName', 'Minion Pro');
exportFigure(fh, struct('SavePath', fullfile('figures', basename), 'Format', 'jpg', 'PaperPosition', [0 0 12, 14]))