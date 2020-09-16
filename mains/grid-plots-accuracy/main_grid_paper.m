%% -- Include stepper files --------------------------------------------------------------------------------------------
addpath(genpath(fullfile(pwd(), '../../stepper/common')));
addpath(fullfile(pwd(), '../../stepper/rk'));
addpath(fullfile(pwd(), '../../stepper/parareal'));
addpath(fullfile(pwd(), 'plotters'));

%% -- Load Equation ----------------------------------------------------------------------------------------------------
equation_name = 'grid-dahlquist';

    cwd = pwd();
    cd(['../../stepper/equations/',equation_name]);
    init
    cd(cwd)

fig_type = 'error';%, 'error' (grid plot showing error of each iteration)  or 'c-region' (shows regios where parareal has converged to fine solution);    
Nbs    = [1 2 4 8]; % number of parareal blocks to use
params = struct(                            ...
    'Nt',           (2^11),                 ... % total steps for experiment
    'num_iters',    0:4,                    ...
    'num_fines',    2.^(2:5),               ...
    'LF',           LF,                     ...
    'NF',           @(t,y) NF(t, y, pars),  ...
    'y0',           y0,                     ...
    'exact',        exact,                  ...
    'z1',           pars.z1,                ...
    'z2',           pars.z2,                ...
    'caxis',        log([1e-16 1e2]),       ...
    'show_speedup_efficiency_title', true,  ...
    'mode',         'exact'                 ...     % "exact" (compares parareal to exact) or "fine" (compares parareal to fine solution)
); 

i_class = 'IMEX';
    
switch(i_class)
    case 'IMEX'
        pr_options = struct(                        ...
            'coarse',       IMEXRKWorker(@IMRK3),   ...
            'fine',         IMEXRKWorker(@IMRK4),   ...
            'num_coarse',   1                       ...
        );
end

params.pr_options = pr_options;

% == Generate Figures ==================================================================================================
for i = 1 : length(Nbs)
    params.Nb = Nbs(i);
    basename  = [i_class, '_', num2str(length(params.z1)), '_Nt_', num2str(params.Nt / params.Nb), '_Nb_', num2str(params.Nb)];
    mat_file  = fullfile('data', ['grid-dahlquist_', basename, '.mat']);
    
    % -- look for data files -------------------------------------------------------------------------------------------
    if(isfile(mat_file))
        load(mat_file, 'y_fine', 'y_parareal', 'y_exact', 'se_parareal');
    else
        [y_fine, y_parareal, y_exact, se_parareal] = ExpGridPlotData(params);
        save(mat_file, 'y_fine', 'y_parareal', 'y_exact', 'se_parareal');
    end
    
    % -- generate and save figure --------------------------------------------------------------------------------------
    switch fig_type
        case 'error'    
            fig_file  = fullfile('figures', ['grid-dahlquist_', fig_type, '-', params.mode, '_', basename]);
            fh = ExpGridPlotPlotter(params, y_fine, y_parareal, y_exact, se_parareal);
            fig_width = 7 * length(params.num_fines);
            fig_height = 7 * length(params.num_iters);
        case 'c-region'
            fig_file  = fullfile('figures', ['grid-dahlquist_', fig_type, '_', basename]);
            fh = DahlquistConvergenceRegionPlotter(params, y_fine, y_parareal, y_exact, se_parareal);
            fig_width = 7 * length(params.num_fines);
            fig_height = 7;
    end
    
    exportFigure(fh, struct('SavePath', fig_file, 'Format', 'jpg', 'PaperPosition', [0 0 fig_width, fig_height]))
    
end