% temporal discretization
tspan = [0 4096];

% parameters 
ng     = 201;
z1_max = .3;
z2_max = .1;
pars   = struct('n', ng, 'z1', linspace(0, z1_max, ng), 'z2', linspace(-z2_max, z2_max, ng), 'epsilon', 0); % 'epsilon', cos(pi/2 - 1/500));

% initial conditions
y0 = ones(1, ng * ng);

% set L and N
LF = L(pars);
NF = @N;

% -- exact solution
[X,Y] = meshgrid(1i * pars.z1, 1i * pars.z2);
exact = @(t) exp((X(:) + Y(:)) * t);