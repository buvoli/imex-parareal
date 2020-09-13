function [R] = rIMRK4D(z1, z2_vec)
%rIMRK4D stability function for LIRK4 method from Eqn (51) from M.P. Calvo, J. de Frutos, J. Novo, "Linearly implicit
%Runge?Kutta methods for advection?reaction?diffusion equations," 2001

% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
% RETURNS
%   R    - value of stability function

[A_im, b_im, A_ex, b_ex] = initIMEXCoefficients();
R = rIMRK(z1, z2_vec, A_im, b_im, A_ex, b_ex);
end

function [A_im, b_im, A_ex, b_ex] = initIMEXCoefficients()
%INITETDRK4 Initializes coefficients for ERK1 for scaler Lambda

A_im = zeros(6);
A_im(2,2) =     1 / 4;
A_im(3,2) =     1 / 2;
A_im(3,3) =     1 / 4;
A_im(4,2) =    17 / 50;
A_im(4,3) =   - 1 / 25;
A_im(4,4) =     1 / 4;
A_im(5,2) =   371 / 1360;
A_im(5,3) = - 137 / 2720; 
A_im(5,4) =   15  / 544;
A_im(5,5) =   1   / 4;     
A_im(6,2) =   25  / 24;
A_im(6,3) = - 49  / 48;
A_im(6,4) =   125 / 16;
A_im(6,5) = - 85  / 12;
A_im(6,6) =   1   / 4;
      
b_im = [0, 25 / 24, - 49 / 48, 125 / 16, - 85 / 12, 1 / 4 ];  

A_ex = zeros(6);
A_ex(2,1) =    1 / 4;
A_ex(3,1) =  - 1 / 4;
A_ex(3,2) =    1;
A_ex(4,1) = - 13 / 100;
A_ex(4,2) =   43 / 75;
A_ex(4,3) =    8 / 75;
A_ex(5,1) =  - 6 / 85;
A_ex(5,2) =   42 / 85;
A_ex(5,3) =  179 / 1360;
A_ex(5,4) = - 15 / 272;
A_ex(6,1) =    0;
A_ex(6,2) =   79 / 24;
A_ex(6,3) =  - 5 / 8;
A_ex(6,4) =   25 / 2;
A_ex(6,5) = - 85 / 6;

b_ex = b_im;

end