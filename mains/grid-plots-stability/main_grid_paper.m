addpath(genpath(fullfile(pwd(), '../../stability'))); % -- include stability folder ------------------------------------
addpath(fullfile(pwd(), 'helpers')); % -- include helper functions -----------------------------------------------------

pars = struct(                                      ...
    'figure_index', 100,                            ...
    'show_unstable_parareal_convergence', true,     ...
    'show_speedup_efficiency_title',      true,     ...
    'scaling_factor', 'total-fine-steps',           ...
    'convergence_norm', Inf,                        ...
    'nfs', 2.^(2:5),                                ...
    'iters', 0:4,                                   ...
    'generate_figures', true,                       ...
    'num_grid_points', 400                          ...
);

%% -- imex integrators -------------------------------------------------------------------------------------------------
domain = struct(                                    ...
    'zoom', {{'very-near'}},                        ...
    'angles', {{'dispersive'}},                     ...
    'ci', {{@rIMRK3}},                              ...
    'fi', {{@rIMRK4}},                              ...
    'nc', 1,                                        ...
    'nt', [512 2048]                                ...
);
multiloop(@MLF, domain, pars);