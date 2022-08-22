classdef (Abstract) FeasibilityDrivenBase < handle
    
   methods (Abstract, Access = public)
       [u, ftstp] = solve(state)
       result = getFeasibilityRegion(obj)
       computeFeasibilityRegion(state, input);
   end 

end