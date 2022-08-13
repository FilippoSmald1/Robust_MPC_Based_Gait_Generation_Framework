classdef (Abstract) FeasibilityDrivenBase
    
   methods (Abstract, Access = public)
       update(obj)
   end 

   properties (Abstract, Access = protected)
       test
   end 
    
end