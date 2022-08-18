classdef RecoveryMode < FeasibilityDrivenBase & handle
    
    methods (Access = public)
        function obj = RecoveryMode(obj)
            obj.a = 2;
            obj.test = 2;
        end
        function [u, ftstp] = solve(obj, state, input)
            u = 0;
            ftstp = 0;
        end
        function result = getFeasibilityRegion(obj)
            result = 0;
        end
        function obj = computeFeasibilityRegion(obj, state, input)
            result = 0;
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