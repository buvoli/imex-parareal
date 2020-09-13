z1_im = linspace(0, 1, 400);
z2_im = linspace(-0.2, 0.2, 400);
convergence_norm = Inf;

pr_params = struct(                             ...
    'np',           64,                         ...
    'nc',           1,                          ...
    'nf',           8,                          ...
    'coarse',       @rIMRK3,                     ...
    'fine',         @rIMRK4,                     ...
    'iterations',   2                           ...
);

%% -- generate data ----------------------------------------------------------------------------------------------------
sf = pr_params.np * pr_params.nf; % set scaling factor based on computational work
amp = @(z1, z2) abs(rPR(z1 * sf, z2 * sf, pr_params));
data_raw_stab = TDStabilityRegionData(amp, z1_im, z2_im, 1i, 1i);       
sv  = @(z1, z2) abs(cnPR(z1 * sf, z2 * sf, pr_params, convergence_norm));
data_raw_svd = TDStabilityRegionData(sv, z1_im, z2_im, 1i, 1i);

%% -- generate plots ---------------------------------------------------------------------------------------------------
data_plot_svd = data_raw_svd;
data_plot_svd(data_raw_svd > 1) = NaN;

data_plot_stab = data_raw_stab;
data_plot_stab(data_raw_stab > 1) = NaN;

data = NaN * ones(size(data_raw_stab));
data(data_raw_stab <= 1 & data_raw_svd <= 1) = 0; % stable method & rapidly convering parareal 
data(data_raw_stab < 1 & data_raw_svd > 1) = 1; % stable method & slowly convering parareal
data(data_raw_stab > 1 & data_raw_svd < 1) = 2; % unstable method & rapidly convering parareal

fh = figure(50);
    surf(z1_im, z2_im, data_plot_svd); 
    view([0 90]); 
    axis([z1_im([1 end]), z2_im([1 end])]);
    shading interp;
    grid off;
    
    sp_p = get(gca,'Position');
    colorbar('Position', [sp_p(1)+sp_p(3)-.065  sp_p(2)+.04  sp_p(3)/12  0.95 * sp_p(4)]);
    
    xlabel('$z_1$','interpreter', 'latex');
    ylabel('$z_2$','interpreter', 'latex');
    title('Convergence Rate -- $\|\mathbf{E}\|$', 'interpreter', 'latex');
    set(gca, 'FontName', 'Minion Pro');
exportFigure(fh, struct('SavePath', 'figures/single-plots/convergence', 'Format', 'pdf', 'PaperPosition', [0 0 1 1]*8))   
     
fh = figure(51);
    contourf(z1_im, z2_im, data_plot_stab, [0 1]); hold on;
    contour(z1_im, z2_im, data_raw_stab, [1 1], 'Color', 'black', 'linewidth', 1);
    hold off;
    colormap([.8 .8 .8; 0 1 0])
    axis([z1_im([1 end]), z2_im([1 end])]);
    xlabel('$z_1$', 'interpreter', 'latex');
    ylabel('$z_2$', 'interpreter', 'latex');
    set(gca, 'FontName', 'Minion Pro');
    title('Stability Region', 'interpreter', 'latex');
exportFigure(fh, struct('SavePath', 'figures/single-plots/stability', 'Format', 'pdf', 'PaperPosition', [0 0 1 1]*8))      
  
fh = figure(52);
    cm = [57 80 151; 80 80 80; 84, 127, 255] / 255;
    surf(z1_im, z2_im, data); shading interp; hold on;
    axis([min(z1_im), max(z1_im), min(z2_im), max(z2_im)]);
    view([0 90])
    colormap(cm)
    grid off 
    contour3(z1_im, z2_im, data_raw_svd * 2, [0 1] * 2, 'Color',  [255 186 65]/255, 'LineWidth', 1);  
    hold off;
    xlabel('$z_1$','interpreter', 'latex');
    ylabel('$z_2$','interpreter', 'latex');
    set(gca, 'FontName', 'Minion Pro');
    title('Stability Convergence Overlay', 'interpreter', 'latex');
 exportFigure(fh, struct('SavePath', 'figures/single-plots/overlay', 'Format', 'pdf', 'PaperPosition', [0 0 1 1]*8))