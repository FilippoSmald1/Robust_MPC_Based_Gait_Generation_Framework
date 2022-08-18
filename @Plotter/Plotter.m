classdef Plotter < handle
    
    methods (Access = public)
        
        function obj = Plotter(logs, input)
            % the constructor is just a setter
            obj.logs = logs;
            obj.input = input;
            obj.rectangle = [obj.input.scheme_parameters.d_z / 2, ...
                             obj.input.scheme_parameters.d_z / 2, ...
                             -obj.input.scheme_parameters.d_z / 2, ...
                             -obj.input.scheme_parameters.d_z / 2, ...
                             obj.input.scheme_parameters.d_z / 2; ...
                             obj.input.scheme_parameters.d_z / 2, ...
                             -obj.input.scheme_parameters.d_z / 2, ...
                             -obj.input.scheme_parameters.d_z / 2, ...
                             obj.input.scheme_parameters.d_z / 2, ...
                             obj.input.scheme_parameters.d_z / 2];
        end
        
        function obj = plotLogs(obj, logs, state)
            
            obj.logs = logs;
            obj.figure_handle = figure(1);
            hold on;
            grid on;
            plot(logs.x_store(1, 1:state.world_time_iter - 1), logs.y_store(1, 1:state.world_time_iter - 1), 'r', 'Linewidth', 2);
            plot(logs.x_store(3, 1:state.world_time_iter - 1), logs.y_store(3, 1:state.world_time_iter - 1), 'b', 'Linewidth', 2);
            legend('CoM', 'ZMP');
            for i = 1 : state.footstep_counter
                c = plot(obj.rectangle(1,:) + logs.actual_footsteps(1,i), obj.rectangle(2,:) + logs.actual_footsteps(2,i), 'm','Linewidth',2, 'Handlevisibility', 'off');
            end
            pbaspect([2 1 1]);
            
        end
        
    end
    
    properties (Access = private)
    
        logs;
        input;
        figure_handle;
        rectangle;
        
    end
    
    
end
    
    