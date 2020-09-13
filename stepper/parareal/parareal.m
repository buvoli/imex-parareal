classdef parareal < handle
    %PARAREAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_coarse
        num_fine
        num_proc
        num_iter
    end
    
    properties (SetAccess = protected)
        coarse_worker
        fine_worker
        h
    end
    
    methods
        
        function this = parareal(options)
            %PARAREAL Construct an instance of this class
            %   pr_params = struct('coarse', CI, 'fine', FI, 'np', np, 'nc', nc, 'nf', nf, 'iterations', iter);
            
            this.coarse_worker = options.coarse;
            this.fine_worker   = options.fine;
            this.num_coarse    = options.num_coarse;
            this.num_fine      = options.num_fine;
            this.num_iter      = options.num_iter;
            this.num_proc      = options.num_proc;
        end
        
        function [t1, y1] = solve(this, L, N, tspan, y0, Nb)
            % SOLVE integrate using multiple parareal blocks
            % L - (array)  \Lambda values L = [\Lambda_1, \Lambda_2, ..., \Lambda_N]
            % N - (handle) nonlinear function
            % tspan - time interval
            % y0 - initial condition
            % Nb - num blocks
        
            h = (tspan(end) - tspan(1)) / (Nb * this.num_proc);
            this.setL(h, L);
            this.setN(N);
            
            t1 = tspan(end);
            y1 = y0;
            
            for i = 1 : Nb
                [t1, y1] = this.step(t1, y1);
            end
            
        end
                 
        function [t1, y1] = step(this, t0, y0)
            %METHOD1 completes one full step of the parareal algorithm
            
            U_k = zeros(length(y0), this.num_proc + 1);
            U_k(:, 1) = y0;
            
            for j = 2 : this.num_proc + 1 % initial coarse guess
                U_k(:, j) = this.coarse_worker.step(t0 + this.h * (j - 1), U_k(:, j-1), this.num_coarse);
            end
            
            U_kp1 = U_k;

            for i = 1 : this.num_iter
                for j = 2 : this.num_proc + 1
                    U_kp1(:, j) = this.coarse_worker.step(t0 + this.h * (j-1), U_kp1(:, j-1), this.num_coarse) + ...
                        this.fine_worker.step(t0 + this.h * (j-1), U_k(:, j-1), this.num_fine) - ...
                        this.coarse_worker.step(t0 + this.h * (j-1), U_k(:, j-1), this.num_coarse);
                end
                U_k = U_kp1;
            end

            t1 = t0 + this.h * this.num_proc;
            y1 = U_kp1(:, end);

        end
        
        function setL(this, h, L)
            %   L   - (array)  of \Lambda values L = [\Lambda_1, \Lambda_2, ..., \Lambda_N]
            this.coarse_worker.setL(h / this.num_coarse, L);
            this.fine_worker.setL(h / this.num_fine, L);
            this.h = h;
        end
        
        function setN(this, N)
            this.coarse_worker.setN(N);
            this.fine_worker.setN(N);
        end
        
    end
end

