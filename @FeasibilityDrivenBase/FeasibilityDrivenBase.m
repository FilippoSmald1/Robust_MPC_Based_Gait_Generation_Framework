classdef (Abstract) FeasibilityDrivenBase
    
   methods (Abstract, Access = public)
       solve(state)
   end 

   properties (Abstract, Access = protected)
       test
   end 
    
end