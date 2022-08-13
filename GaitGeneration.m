classdef GaitGeneration < FeasibilityDrivenBase
    
    methods (Access = public)
        function obj = GaitGeneration(obj)
            obj.a = 2;
            obj.test = 2;
        end
        function obj = update(obj)
            
        end
    end
    
    %properties (Access = public)
    %    a = 1;
    %end
    
    properties (Access = private)
        a
    end
    
    properties (Access = protected)
        test = 1;
    end
    
end