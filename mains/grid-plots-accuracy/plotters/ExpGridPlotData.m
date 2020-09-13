function [y_fine, y_parareal, y_exact, se_parareal] = ExpGridPlotData(params)
%UNTITLED Summary of this function goes here
%   struct with fields
%
%   pr_params   (struct) - parareal parameters
%   LF          () - linear operator
%   NF          handle - nonlinear function

num_fines = params.num_fines;   % num procs
num_iters = params.num_iters;   % num iterations
Nt        = params.Nt;          % total number of timesteps
Nb        = params.Nb;          % number of blocks
tspan     = [0 Nt];

% == compute fine solution =============================================================================================
pr_options = params.pr_options;
pr_options.num_iter = 1;
pr_options.num_proc = 1;
pr_options.num_fine = Nt;
PR = parareal(pr_options);
[~, y_fine] = PR.solve(params.LF, params.NF, tspan, params.y0, 1);

y_exact = params.exact(tspan(end));

% == compute parareal solution =========================================================================================
y_parareal = cell(length(num_iters), length(num_fines));
se_parareal = cell(length(num_iters), length(num_fines));
for i = 1 : length(num_iters)
    for j = 1 : length(num_fines)
        
        % 1. check that Np divides Nt
        if(mod(Nt, Nb) ~= 0)
            error('Nb must divide Nt')
        end

        % 2. check that NF divides Nt
        if(mod(Nt/Nb, num_fines(j)) ~= 0)
            error('num_fines must divide Nt')
        end
        
        % 2. set up parareal params
        pr_options = params.pr_options;
        pr_options.num_iter = num_iters(i);
        pr_options.num_proc = (Nt / Nb) / num_fines(j); % processors per block
        pr_options.num_fine = num_fines(j);
        
        % 3. run parareal
        PR = parareal(pr_options);
        [~, y_parareal{i,j}] = PR.solve(params.LF, params.NF, tspan, params.y0, Nb);
        
        % 4. compute speedup and efficiency of method
        [s, e] = pararealSpeedup(pr_options);
        se_parareal{i,j} = [s, e];
        
    end
end

end