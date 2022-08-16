classdef RobustGaitGenerationScheme

    methods (Access = public)
        function obj = RobustGaitGenerationScheme(obj)
            obj.a = 2;
            obj.test = 2;
            obj.sm_instance = StandardMode;
            obj.rm_instance = RecoveryMode;
        end
        function obj = update(obj)
            
        end
    end
    
    properties (Access = public)
        
        a = 1;
        test = 1;
        
    end

    properties (Access = private)
        
        sm_instance;
        rm_instance;
        
    end
    
end