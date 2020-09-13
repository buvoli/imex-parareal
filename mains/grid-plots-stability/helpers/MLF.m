function MLF(parameters)
%MLF Multiloop Function - function that should be passed to multiloop to generate figures.

% extract function names 
foo = functions(parameters.ci);
ci_name = foo.function;
foo = functions(parameters.fi);
fi_name = foo.function;

% exit if course integrator has higher order then fine 
if(ci_name(end) > fi_name(end))
    return
end

disp(parameters);

% determine integrator type
if(contains(ci_name, 'ERK'))
    i_class = 'exp';
elseif(contains(ci_name, 'IMRK'))
    i_class = 'imex';
else
    i_class = 'unknown';
end

num_grid_points = parameters.num_grid_points;
[z1_r, z2_r, z1_angle, z2_angle] = z1z2(parameters.zoom, parameters.angles, num_grid_points);
parameters.z1_r = z1_r;
parameters.z2_r = z2_r;
parameters.z1_angle = z1_angle;
parameters.z2_angle = z2_angle;

% -- load data if existant ---------------------------------------------------------------------------------------------
parent_dir = fullfile('figures', ['nt-', num2str(parameters.nt)]); % organize by num fine
ckmkdir(parent_dir);
filename = fullfile(parent_dir, join({parameters.angles, parameters.zoom, num2str(num_grid_points),  i_class, 'course', ci_name, 'nc', num2str(parameters.nc), 'fine', fi_name, 'nt', num2str(parameters.nt)}, '_'));
data_file = [filename{:}, '.mat'];

if(isfile(data_file))
    load(data_file, 'data_stab', 'data_svd', 'data_se');
else
    [data_stab, data_svd, data_se] = pararealGridData(parameters); 
    save(data_file, 'parameters', 'data_stab', 'data_svd', 'data_se')
end

if(parameters.generate_figures)
    fh = pararealGridPlotter(data_stab, data_svd, data_se, parameters);
    fig_width = 7 * length(parameters.nfs);
    fig_height = 7 * length(parameters.iters);
    exportFigure(fh, struct('SavePath', filename, 'Format', 'jpg', 'PaperPosition', [0 0 fig_width, fig_height]))
end

end

