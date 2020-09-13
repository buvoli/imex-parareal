function [s] = cnPR(z1, z2, params, convergence_norm)
%CNPR returns the norm of the RK-parareal iteration
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
%   params (struct) - parameters used for parareal. Must contain fields
%       nc (integer) - number of coarse steps
%       nf (integer) - number of fine steps per coarse step
%       np (integer) - number of processors
%       coarse (@(z1, z2)) - coarse integrator stability function
%       fine (@(z1, z2)) - fine integrator stability function
%   convergence_norm (integer - optional) - norm used for determining convergence region.
% RETURNS
%   s - largest singular value of iteration matrix

if(nargin < 4)
    convergence_norm = 2;
end

nc = params.nc; % number of total steps
np = params.np; % number of coarse steps per processor
nf = params.nf; % number of fine steps per processor

C = params.coarse(z1 / (nc * np), z2 / (np * nc)) .^ (nc);
F = params.fine(z1 / (np * nf), z2 / (np * nf)).^(nf);
I = eye(np+1);
    
s = zeros(size(z2));
for i = 1 : length(z2)
    
    switch(convergence_norm)
        case Inf % use explicit form for inf norm
            if(abs(C(i)) ~= 1)
                s(i) = ((1 - abs(C(i))^(np)) / (1 - abs(C(i)))) * abs(C(i) - F(i));
            else
                s(i) = np * abs(C(i) - F(i));
            end
        otherwise % generate iteration matrix and compute norm
            Mc = diag(ones(np+1,1),0) - C(i) * diag(ones(np, 1), -1);
            Mf = diag(ones(np+1,1),0) - F(i) * diag(ones(np, 1), -1);
            PM = I - Mc \ Mf;
            if(any(isnan(PM), 'all'))
                s(i) = NaN;
            else
                s(i) = norm(PM, convergence_norm);
            end        
    end
    
end
end