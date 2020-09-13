classdef IMEXRKWorker < AbstractRKWorker
    %RKWORKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        % rk matrices
        A_im
        b_im
        A_ex
        b_ex
        c
        nstages
        % method
        coeffGenerator
    end
    
    properties (Access = private)
        Y          % stages
        YN         % nonlinear stage derivatives
        YL         % linear stage derivatives
    end
    
    methods
        
        function this = IMEXRKWorker(coeffGenerator)
            %RKWORKER Construct an instance of this class
            %   coeffGenerator (handle) function that returns the method coefficients
            [this.A_im, this.b_im, this.A_ex, this.b_ex, this.c] = coeffGenerator();
            this.validateCoefficients()
            this.nstages = size(this.A_im, 2);
            this.coeffGenerator = coeffGenerator;
        end
        
        function [y0, t0] = step(this, t0, y0, Nt)
            %METHOD1 computes Nt steps of the RK method
            
            % read variables into local variables for faster speed
            Y = this.Y;
            YL = this.YL;
            YN = this.YN;
            A_im = this.A_im;
            b_im = this.b_im;
            A_ex = this.A_ex;
            b_ex = this.b_ex;            
            c = this.c;
            nstages = this.nstages;
            h = this.h;
            L = this.L;
            
            for n = 1 : Nt
                % -- compute stage values ------------------------------------------------------------------------------------------
                for i = 1 : nstages
                    Y(:, i) = y0;
                    for j = 1 : (i - 1)
                        Y(:, i) = Y(:, i) + (h * A_im(i,j)) * (YL(:, j)) + (h * A_ex(i,j)) * YN(:, j);
                    end
                    Y(:, i) = Y(:, i) ./ ( 1 - (h * A_im(i,i)) * L);
                    % -- evaluate stage derivatives
                    YL(:,i) = L .* Y(:, i);
                    YN(:,i) = this.N(t0 + h * c(i), Y(:,i));
                end
                % -- compute output ------------------------------------------------------------------------------------------------
                for i = 1 : nstages
                    y0 = y0 + (h * b_im(i)) .* YL(:, i) + (h * b_ex(i)) .* YN(:, i);
                end
                t0 = t0 + h;   
            end
            
        end
        
        function setL(this, h, L)
            %   h (double) stepsize
            %   L (array)  of \Lambda values L = [\Lambda_1, \Lambda_2, ..., \Lambda_N]
     
            this.L = L(:);
            this.h = h;
            this.Y  = zeros(length(L), this.nstages); 
            this.YN = zeros(length(L), this.nstages);
            this.YL = zeros(length(L), this.nstages);
        end
        
        function setN(this, N)
            % N (handle) nonlinear function
            
            this.N = N;
        end
        
    end
    
    methods (Access = protected)
        
        function validateCoefficients(this)

        s = size(this.A_im, 1); % number of total stages (including explicit stages)
        dims = [size(this.A_im, 2) size(this.A_ex, 1) size(this.A_ex, 2) length(this.b_im) length(this.b_ex) length(this.c)];

        if(~all(dims == s))
            error('invalid coefficient matrices: A_im, A_ex should be sxs while b_im, b_ex, and c should have length s')
        end

        end
        
    end
end

