function [ S, E ] = pararealSpeedup(options)
%SPEEDUP Returns speedup and efficiencry for parareal method
%   Detailed explanation goes here

cost_coarse = cost(options.coarse);
cost_fine   = cost(options.fine);
    
alpha = (cost_coarse * options.nc) / (cost_fine * options.nf);
K     = options.iterations;
Np    = options.np;

S = Np / ( Np * alpha + K * (1 + alpha));
E = S / options.np;

end

function c = cost(handle)
    handle_info = functions(handle);
    handle_name = handle_info.function;
    switch(handle_name)
        case {"rERK1", "rIMRK1"}
            c = 1;
        case {"rERK2", "rIMRK2"}
            c = 2;
        case {"rERK3", "rIMRK3"}
            c = 3;
        case {"rERK4", "rIMRK4"}
            c = 4;
        otherwise
            warning('Unknown integator; Cannot determine cost.')
            c = NaN;
    end
end

