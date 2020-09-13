function [r] = rPR(z1, z2, params)
%RPR returns the stability function for parareal run using a course and fine RK method
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
%   params (struct) - parameters used for parareal. Must contain fields
%       nc (integer) - number of coarse steps
%       nf (integer) - number of fine steps per coarse step
%       iterations (integer) - number of parareal iterations
% RETURNS
%   r - value of stability function

nc = params.nc; % number of total steps
np = params.np; % number of processors or blocks
nf = params.nf; % number of fine steps per coarse step
it = params.iterations;

C = params.course(z1 / (nc * np), z2 / (np * nc)) .^ (nc);
F = params.fine(z1 / (np * nf), z2 / (np * nf)) .^ (nf);

c1 = [1; zeros(np, 1)];
c2 = [zeros(1, np), 1];
D  = ones(np + 1, 2);

r = zeros(size(z2));
for i = 1 : length(z2)
    
    % 1. generate sparse Mg Matrix -------------------------------------------------------------------------------------
    D(:, 1) = C(i);
    Mg = spdiags(D, [-1, 0], np + 1, np + 1); 
    % 2. generate sparse Mf Matrix -------------------------------------------------------------------------------------
    D(:, 1) = F(i);
    Mf = spdiags(D, [-1, 0], np + 1, np + 1);    
    % 3. compute c_2 ( \sum_{j=0}^k E^j ) M_g^{-1} c_1 where E = I - Mg^{-1} * Mf --------------------------------------
    x0 = Mg \ c1;
    r(i) = c2 * x0;                     % zeroth term in sum
    for j = 1 : min(it, np * nf)
        x0 = x0 - (Mg \ (Mf * x0));     % compute ( E^j ) M_g^{-1} c_1
        r(i) = r(i) + c2 * x0;          % add on c2 * ( E^j ) M_g^{-1} c_1
    end
    
end
end