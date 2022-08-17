classdef BuildRestrictionFunction
   
    methods (Access = public)
        
        function obj = BuildRestrictionFunction(input)
            
            obj.input = input;
            costFunction = @(x) x(1)^2 + (x(2) - obj.input.scheme_parameters.T_c)^2;            
            x0 = [obj.input.scheme_parameters.d_z / (2.0 * obj.input.scheme_parameters.T_c), ...
                  obj.input.scheme_parameters.T_c];
              
            obj.input.scheme_parameters.gamma_max = (2 * obj.input.scheme_parameters.dist_range(1,1) / obj.input.scheme_parameters.eta ^ 2) * ...
                                                    (exp(obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) - 1) ...
                                                    * (1 + obj.input.scheme_parameters.alpha);
            nlcon = @(x) constraintFunction (obj, x);                               
            x = fmincon(costFunction, x0, [], [], [], [], [], [], nlcon);
            obj.restriction_x = zeros(obj.input.scheme_parameters.T_c / obj.input.scheme_parameters.delta, 1);
            for  i = 1 : (obj.input.scheme_parameters.T_c / obj.input.scheme_parameters.delta) - 1
                if i * obj.input.scheme_parameters.delta < x(2)
                    obj.restriction_x(i+1) = x(1) * (i * obj.input.scheme_parameters.delta - obj.input.scheme_parameters.delta);
                else
                    obj.restriction_x(i+1) = x(1) * (x(2) - obj.input.scheme_parameters.delta);
                end
            end
            
            obj.input.scheme_parameters.gamma_max = (2 * obj.input.scheme_parameters.dist_range(2,1) / obj.input.scheme_parameters.eta ^ 2) * ...
                                                    (exp(obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) - 1) ...
                                                    * (1 + obj.input.scheme_parameters.alpha);
            nlcon = @(x) constraintFunction (obj, x);                                  
            x = fmincon(costFunction, x0, [], [], [], [], [], [], nlcon);
            obj.restriction_y = zeros(obj.input.scheme_parameters.T_c / obj.input.scheme_parameters.delta, 1);
            for  i = 1 : (obj.input.scheme_parameters.T_c / obj.input.scheme_parameters.delta) - 1
                if i * obj.input.scheme_parameters.delta < x(2)
                    obj.restriction_y(i+1) = x(1) * (i * obj.input.scheme_parameters.delta - obj.input.scheme_parameters.delta);
                else
                    obj.restriction_y(i+1) = x(1) * (x(2) - obj.input.scheme_parameters.delta);
                end
            end
            
        end
          
        function restriction_x = getRestrictionX(obj)
            restriction_x = obj.restriction_x;
        end
        
        function restriction_y = getRestrictionY(obj)
            restriction_y = obj.restriction_y;
        end
        
        function obj = plotRestrictions(obj)
            figure
            clf;
            hold on;
            grid on;
            plot(obj.restriction_x,'Linewidth',2);
            plot(obj.restriction_y,'Linewidth',2);
            legend('res x', 'res y');
            pbaspect([2 1 1]);
        end
        
    end
    
    methods (Access = private)
        
        function [c,ceq] = constraintFunction(obj, x)

            c =  [-x(1); ...
                 x(1)*(x(2) - obj.input.scheme_parameters.delta) - obj.input.scheme_parameters.d_z / 2.0 + obj.input.scheme_parameters.epsilon];
            ceq = [(x(1) / obj.input.scheme_parameters.eta) * (exp(- obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) - 1) * ...
                   ((obj.input.scheme_parameters.eta * (obj.input.scheme_parameters.delta - x(2) -1))*exp(- obj.input.scheme_parameters.eta*x(2)) - 1)  + ...
                   x(1) * (x(2) - obj.input.scheme_parameters.delta) * (exp(obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) - 1) * ...
                   (exp(-obj.input.scheme_parameters.eta * x(2)) - exp(-obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_c)) - ...
                   obj.input.scheme_parameters.mi_max - obj.input.scheme_parameters.gamma_max];
           
         end     
    end
    
    properties (Access = private)
        
        input;
        restriction_x;
        restriction_y;
        
    end
    
end