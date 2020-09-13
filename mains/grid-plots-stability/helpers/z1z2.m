function [z1_r, z2_r, z1_angle, z2_angle] = z1z2(zoom_str, angle_str, np)
% Z1Z2 produces imaginary ranges for z1 and z2

if(nargin == 1)
    np = 200;
end

% z1_L is length of z1 interval [0, z1_L]
% z2_hL is half-length of z2 interval [-z2_l, z2_l]
switch(zoom_str)
    case 'very-near'
        z2_hL = .1;
        z1_L = .3;
    case 'far'
        z2_hL = 4;
        z1_L = 60;
    otherwise
        z2_hL = 1.6;
        z1_L = 16;
end

z1_r = linspace(0, z1_L, np);
z2_r = linspace(-z2_hL, z2_hL, np);

switch(angle_str)
    case 'dispersive'
        z1_angle = 1i;
        z2_angle = 1i;
    case 'damped-dispersive'
        a_eps = 1 / 500;
        z1_angle = exp(1i * (pi/2 + a_eps));
        z2_angle = exp(1i * (pi/2 - a_eps));
    case 'dissipative'
        z1_angle = -1;
        z2_angle = -1;
end

end