classdef AbstractRKWorker < handle
    %RKWORKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        % problem properties
        L
        N
        % other properties
        h
    end
        
    methods (Abstract)
                
        step(this, t0, y0, Nt)
        
        setL(this, h, L)
        
        setN(this, N)
        
    end
            
end

