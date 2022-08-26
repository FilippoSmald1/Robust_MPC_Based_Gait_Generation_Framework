classdef Plotter < handle
    
    methods (Access = public)
        
        function obj = Plotter(logs, input)
            % the constructor is just a setter
            obj.logs = logs;
            obj.input = input;
            obj.color = [0.4660 0.6740 0.1880];
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
            clf;
            hold on;
            grid on;
            com = plot(logs.x_store(1, 1:state.world_time_iter - 1), logs.y_store(1, 1:state.world_time_iter - 1), 'r', 'Linewidth', 2);
            zmp = plot(logs.x_store(3, 1:state.world_time_iter - 1), logs.y_store(3, 1:state.world_time_iter - 1), 'b', 'Linewidth', 2);
            for i = 1 : state.footstep_counter
                c = plot(obj.rectangle(1,:) + logs.actual_footsteps(1,i), obj.rectangle(2,:) + logs.actual_footsteps(2,i), 'm','Linewidth',2, 'Handlevisibility', 'off');
            end
            rectangle_feasibility = [state.feasibility_region(2,1), state.feasibility_region(2,1), state.feasibility_region(1,1), state.feasibility_region(1,1), state.feasibility_region(2,1); ...
                                     state.feasibility_region(4,1), state.feasibility_region(3,1), state.feasibility_region(3,1), state.feasibility_region(4,1), state.feasibility_region(4,1)];

            feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color);
            feas_region_patch.FaceAlpha = 0.3;
            feas_region_patch.EdgeColor = [0.4660 0.6740 0.1880];
            dcm = scatter(state.x(1,1) + state.x(2,1)/ obj.input.scheme_parameters.eta, ...
                        state.y(1,1) + state.y(2,1)/ obj.input.scheme_parameters.eta,...
                        'g','Linewidth',2);
            axis equal
            axis([-0.1 2.7 -0.25 0.55])   
            pbaspect([2 1 1]);
            legend([com, zmp, feas_region_patch, dcm], {'CoM', 'ZMP', 'current feas region', 'current dcm'});
            
        end
        
    end
    
    properties (Access = private)
    
        logs;
        input;
        figure_handle;
        rectangle;
        color;
        
    end
    
    
end
    
    