classdef (Abstract) FeasibilityDrivenBase
    
   methods (Abstract, Access = public)
       solve(obj)
   end 

   properties (Abstract, Access = protected)
       test
   end 
    
end