function [ fh ] = DahlquistConvergenceRegionPlotter(params, y_fine, y_parareal, y_exact, ~)
%UNTITLED Summary of this function goes here
%   params
%       caxis - colorbar axis
%       z1 - axis for z1
%       z2 - axis for z2

% -- numerical parameters ----------------------------------------------------------------------------------------------
mode      = params.mode;
num_fines = params.num_fines;   % num procs
num_iters = params.num_iters;   % num iterations
len_num_iters = length(num_iters);
len_num_fines = length(num_fines);

fh = figure();


error_fine = abs(y_fine - y_exact); % error of fine solution
% -- parareal ----------------------------------------------------------------------------------------------------------
for j = 1 : len_num_fines
    ax = subplot(1, len_num_fines, j);
    cregions = NaN * zeros(size(y_fine));
    
    for i = len_num_iters:-1:1
        error = abs(y_fine - y_parareal{i,j}); % error of kth iteration
        cregions(error <= error_fine + 100 * eps) = num_iters(i) + 1/2; % add eps in case solution is closer than machine precision (or gaps will form since parareal solution will never converge past machine precision until k=np)
    end
    
    convPlot(cregions, params)
    
    title_str = ['{\boldmath$N_F = ', num2str(num_fines(j)), ',~ N_p=', num2str(params.Nt / (params.Nb * num_fines(j))), '$}'];
    title({title_str, ''}, 'interpreter', 'latex');
    set(ax, 'FontName', 'Minion Pro');
    
    if(j == len_num_fines)
        sp_p = get(ax, 'Position');
        colorbar( ...
             'Position', [sp_p(1)+sp_p(3)+0.01  sp_p(2)  sp_p(3)/8  sp_p(4)], ...
             'XTick', num_iters + 1/2, ...
             'XTickLabel', arrayfun(@(x) ['k=', num2str(x)], num_iters, 'UniformOutput', false))
        caxis([min(num_iters) max(num_iters) + 1]);
    end
      
end

end

function convPlot(cregions, params)
    z1 = params.z1;
    z2 = params.z2;
    nz1 = length(params.z1);
    nz2 = length(params.z2);
    num_iters = params.num_iters;   % num iterations
    len_num_iters = length(num_iters);
    surf(params.z1, params.z2, reshape(cregions, nz1, nz2)); 
    shading interp
    axis([min(z1) max(z1) min(z2) max(z2) min(num_iters) max(num_iters) + 1/2])
    colormap(linspace(0, .9, len_num_iters).' .* [255, 255, 255] / 255);
    view([0 90]);
    grid off;
    
end