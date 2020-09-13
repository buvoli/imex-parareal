function [data_raw_stab, data_raw_svd, data_raw_se] = pararealGridData(options)
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

% -- Read Parameters ---------------------------------------------------------------------------------------------------
z1_r = options.z1_r;
z1_angle = options.z1_angle;
z2_r = options.z2_r;
z2_angle = options.z2_angle;

nc    = options.nc;         % num course steps
nt    = options.nt;         % num time steps
nfs   = options.nfs;        % num fine steps
iters = options.iters;      % num iters
ci    = options.ci;
fi    = options.fi;
c_nrm = options.convergence_norm;

% -- generate data -----------------------------------------------------------------------------------------------------
num_nf = numel(nfs);
num_iter = numel(iters);
data_raw_svd  = cell(num_nf, num_iter);
data_raw_stab = cell(num_nf, num_iter);
data_raw_se   = cell(num_nf, num_iter);

for i = 1 : num_nf
    if(mod(nt, nfs(i)) ~= 0)
        error('nf must divide nt')
    end
    np = nt / nfs(i);
    pr_params = struct('course', ci, 'fine', fi, 'np', np, 'nc', nc, 'nf', nfs(i));
    switch(options.scaling_factor)
        case 'total-fine-steps'
            sf = nt; % set scaling factor based on computational work of fine integrator
        otherwise
            error('invalid scaling factor.');
    end
    for j = 1 : num_iter
        pr_params.iterations = iters(j);
        amp = @(z1, z2) abs(rPR(z1 * sf, z2 * sf, pr_params));
        data_raw_stab{i,j} = TDStabilityRegionData(amp, z1_r, z2_r, z1_angle, z2_angle);
        [s, e] = pararealSpeedup(pr_params);
        data_raw_se{i,j} = [s e];
    end
    sv  = @(z1, z2) abs(cnPR(z1 * sf, z2 * sf, pr_params, c_nrm));
    data_raw_svd{i} = TDStabilityRegionData(sv, z1_r, z2_r, z1_angle, z2_angle); 
end

end