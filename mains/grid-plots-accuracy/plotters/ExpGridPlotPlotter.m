function [ fh ] = ExpGridPlotPlotter(params, y_fine, y_parareal, y_exact, se_parareal)
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

% -- fine solution -----------------------------------------------------------------------------------------------------
for j = 1 : len_num_fines
    ax = subplot(len_num_iters + 1, len_num_fines, j);
    errorPlot(y_fine - y_exact, params)
    if(j == 1)
        ylabel({'{\boldmath $k = N_p$}', ''}, 'interpreter', 'latex');
    end
    title_str = ['{\boldmath$N_F = ', num2str(num_fines(j)), ',~ N_p=', num2str(params.Nt / (params.Nb * num_fines(j))), '$}'];
    title({title_str, ''}, 'interpreter', 'latex');
    if(j == len_num_fines)
        sp_p = get(ax,'Position');
        colorbar('Position', [sp_p(1)+sp_p(3)+0.01  sp_p(2)  sp_p(3)/8  sp_p(4)]);
    end
    set(ax, 'FontName', 'Minion Pro');
end

% -- parareal ----------------------------------------------------------------------------------------------------------
for i = 1 : len_num_iters
    for j = 1 : len_num_fines
            ax = subplot(len_num_iters + 1, len_num_fines, len_num_fines * i + j);
            
            switch mode
                case "exact"
                    errorPlot(y_parareal{i,j} - y_exact, params);
                case "fine"
                    errorPlot(y_parareal{i,j} - y_fine, params);
            end
            
            if(j == 1)
                ylabel({['{\boldmath $k = ', num2str(num_iters(i)), '$}'], ''}, 'interpreter', 'latex');
            end
            if(params.show_speedup_efficiency_title)
                title(sprintf('S: %.3f,   E: %.3f', se_parareal{i, j}), 'interpreter', 'latex')
            end
            set(ax, 'FontName', 'Minion Pro');
    end
end

end

function errorPlot(error, params)
    nz1 = length(params.z1);
    nz2 = length(params.z2);
    surf(params.z1, params.z2, log(abs(reshape(error, nz1, nz2)) + eps));
    shading interp
    caxis(params.caxis)
    grid off;
    view([0 90]);
end