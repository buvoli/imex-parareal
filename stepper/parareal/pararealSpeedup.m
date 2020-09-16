function [ S, E ] = pararealSpeedup(options)
%SPEEDUP Returns speedup and efficiencry for parareal method
%   Detailed explanation goes here

cost_coarse = cost(options.coarse);
cost_fine   = cost(options.fine);
    
alpha = (cost_coarse * options.num_coarse) / (cost_fine * options.num_fine);
K     = options.num_iter;
Np    = options.num_proc;

S = Np / ( Np * alpha + K * (1 + alpha));
E = S / options.num_proc;

end

function c = cost(rk_worker)
    handle_info = functions(rk_worker.coeffGenerator);
    handle_name = handle_info.function;
    switch(handle_name)
        case {"ERK1", "IMRK1"}
            c = 1;
        case {"ERK2", "IMRK2"}
            c = 2;
        case {"ERK3", "IMRK3"}
            c = 3;
        case {"ERK4"}
            c = 4;
        case {"IMRK4"}
            c = 5;
        otherwise
            warning('Unknown integator; Cannot determine cost.')
            c = NaN;
    end
end

