classdef lleProgress
    
    properties
        fig
        t_evol
    end
    
    methods
        function obj = lleProgress(t_evol)
            progressbar('solving...');
            obj.t_evol = t_evol;
        end
        
        function report(obj, t, ~)
            progressbar(t/obj.t_evol);
        end
        
    end
    
end