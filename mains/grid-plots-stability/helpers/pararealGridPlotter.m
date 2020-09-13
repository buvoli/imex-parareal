function [fh] = pararealGridPlotter(data_raw_stab, data_raw_svd, data_raw_se, options)
%PARAREALGRID produces a grid of plots that show parareal stability and convergence regions. Number of processors is
%rendered along the x axis and number of iterations along the y axis.
% PARAMETERS
%   options (struct) with the following fields
%       z1_im: vector   - imaginary z1 points
%       z2_im: vector   - imaginary z2 points
%       nc: number      - number of coarse integrator steps
%       nf: number      - number of fine integrator steps
%       ci: @(z1,z2)    - stability function for course integrator
%       fi: @(z1,z2)    - stability function for fine integrator
%       nps: vector     - different values of processors that should be considered.
%       iters: vector   - different values of iterations that should be considered.
%       show_unstable_parareal_convergence: boolean - if true, regions where parareal converges but integrator diverges 
%                                                     will be labaled in red.
%       figure_index: number - 

% -- Read Parameters ---------------------------------------------------------------------------------------------------
z1_r = options.z1_r;
z2_r = options.z2_r;
nfs   = options.nfs;        % num fine steps
iters = options.iters;      % num iters

% -- generate data -----------------------------------------------------------------------------------------------------
num_fs = numel(nfs);
num_iter = numel(iters);

% -- Generate Plots ----------------------------------------------------------------------------------------------------
if(options.show_unstable_parareal_convergence)
    cm = [57 80 151; 80 80 80; 63 111 255] / 255; % colormap
else
    cm = [77 77 77; 153 153 153] / 255; % colormap
end

fh = figure(options.figure_index);
for j = 1 : num_fs
    for i = 1 : num_iter + 1      
        ax = subplot(num_iter + 1, num_fs , (i - 1) * num_fs + j);
        if(i == 1) % show svd amp in first row
            pr_svd = data_raw_svd{j};
            pr_svd(pr_svd > 1) = NaN;            
            surf(z1_r, z2_r, pr_svd); shading interp;
            axis([z1_r(1) z1_r(end) z2_r(1) z2_r(end) 0 1]);
            view([0 90])
            colormap(ax, parula)
            grid off;
            if(j == num_fs)
                sp_p = get(ax,'Position');
                colorbar('Position', [sp_p(1)+sp_p(3)+0.01  sp_p(2)  sp_p(3)/8  sp_p(4)]);
            end
        else
            pr_svd = data_raw_svd{j};
            pr_svd(pr_svd > 1e30) = 1e30; % clip large values to prevent plot error "Limits are too large"
            
            di = i - 1; % shift data index by one
            data = NaN * ones(size(data_raw_stab{j,di}));
            data(data_raw_stab{j,di} <= 1 & data_raw_svd{j} <= 1) = 0; % stable method & rapidly convering parareal 
            data(data_raw_stab{j,di} < 1 & data_raw_svd{j} > 1) = 1; % stable method & slowly convering parareal

            if(options.show_unstable_parareal_convergence)
                data(data_raw_stab{j,di} > 1 & data_raw_svd{j} < 1) = 2; % unstable method & rapidly convering parareal
            end

            surf(z1_r, z2_r, data); shading interp; hold on;
            axis([min(z1_r), max(z1_r), min(z2_r), max(z2_r)]);
            view([0 90])
            colormap(ax, cm)
            grid off 
            contour3(z1_r, z2_r, pr_svd * 2, [0 1] * 2, 'Color',  [255 186 65]/255, 'LineWidth', 1); hold off;
            
            if(options.show_speedup_efficiency_title)
                title(sprintf('S: %.3f,   E: %.3f', data_raw_se{j, di}), 'interpreter', 'latex')
            end
            
        end  
        
        if(i == 1)
            title_str = ['{\boldmath$N_F = ', num2str(nfs(j)), ',~ N_p=', num2str(options.nt / (nfs(j))), '$}'];
            title({title_str, ''}, 'interpreter', 'latex');  
        end
        
        if(j == 1 && i == 1)
            switch options.convergence_norm
                case 2
                    title_str = '{\boldmath $\max_j \sigma_j$}';
                case Inf
                    title_str = '{\boldmath $\|E\|_{\infty}$}';
                otherwise
                    title_str = sprintf('{\\boldmath $\|E\|_{%s}$}', num2str(options.convergence_norm));
            end
            ylabel({title_str, ''}, 'interpreter', 'latex');  
        end
        if(j == 1 && i > 1)
            title_str = sprintf('{\\boldmath $k = %s$}', num2str(iters(i - 1))); 
            ylabel({title_str, ''}, 'interpreter', 'latex');  
        end
        set(ax, 'FontName', 'Minion Pro');
        
    end
    
end

end