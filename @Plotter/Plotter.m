classdef Plotter < handle
    
    methods (Access = public)
        
        function obj = Plotter(logs, input, simulation_parameters)
            % the constructor is just a setter
            obj.logs = logs;
            obj.sim_parameters = simulation_parameters;
            obj.input = input;
            obj.color = [0.4660 0.6740 0.1880];
            obj.color_estimate = [0.2660 0.3740 0.8880];
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
            obj.initial_ds = [obj.input.scheme_parameters.d_z / 2, ...
                             obj.input.scheme_parameters.d_z / 2, ...
                             -obj.input.scheme_parameters.d_z / 2, ...
                             -obj.input.scheme_parameters.d_z / 2, ...
                             obj.input.scheme_parameters.d_z / 2; ...
                             0.09+obj.input.scheme_parameters.d_z / 2, ...
                             -0.09-obj.input.scheme_parameters.d_z / 2, ...
                             -0.09-obj.input.scheme_parameters.d_z / 2, ...
                             0.09+obj.input.scheme_parameters.d_z / 2, ...
                             0.09+obj.input.scheme_parameters.d_z / 2];                         
        end
        
        function obj = plotLogs(obj, logs, state)
            
            obj.counter = obj.counter + 1;
            obj.logs = logs;
            obj.figure_handle = figure(1);
            clf;
            hold on;
            grid on;
            com = plot(logs.x_store(1, 1:state.sim_iter - 1), logs.y_store(1, 1:state.sim_iter - 1), 'r', 'Linewidth', 2);
            zmp = plot(logs.x_store(3, 1:state.sim_iter - 1), logs.y_store(3, 1:state.sim_iter - 1), 'b', 'Linewidth', 2);
            
            c = plot(obj.initial_ds(1,:), obj.initial_ds(2,:), 'm', 'Linewidth',2, 'Handlevisibility', 'off');
            
            if state.footstep_counter >= 1
                for i = 2 : state.footstep_counter
                    c = plot(obj.rectangle(1,:) + logs.actual_footsteps(1,i), obj.rectangle(2,:) + logs.actual_footsteps(2,i), 'm','Linewidth',2, 'Handlevisibility', 'off');
                end
            end
  
            rectangle_feasibility = [state.feasibility_region(2,1), state.feasibility_region(2,1), state.feasibility_region(1,1), state.feasibility_region(1,1), state.feasibility_region(2,1); ...
                                     state.feasibility_region(4,1), state.feasibility_region(3,1), state.feasibility_region(3,1), state.feasibility_region(4,1), state.feasibility_region(4,1)];
                                 
            rectangle_feasibility_estimate = [state.feasibility_region(6,1), state.feasibility_region(6,1), state.feasibility_region(5,1), state.feasibility_region(5,1), state.feasibility_region(6,1); ...
                                     state.feasibility_region(8,1), state.feasibility_region(7,1), state.feasibility_region(7,1), state.feasibility_region(8,1), state.feasibility_region(8,1)];
                                 
            feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color);
            feas_region_patch.FaceAlpha = 0.3;
            feas_region_patch.EdgeColor = [0.4660 0.6740 0.1880];
            
            rectangle_feasibility = [state.feasibility_region(2,2), state.feasibility_region(2,2), state.feasibility_region(1,2), state.feasibility_region(1,2), state.feasibility_region(2,2); ...
                                     state.feasibility_region(4,2), state.feasibility_region(3,2), state.feasibility_region(3,2), state.feasibility_region(4,2), state.feasibility_region(4,2)];
            
            feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color);
            feas_region_patch.FaceAlpha = 0.3;
            feas_region_patch.EdgeColor = [0.2660 0.3740 0.8880];   
            
            rectangle_feasibility = [state.feasibility_region(2,3), state.feasibility_region(2,3), state.feasibility_region(1,3), state.feasibility_region(1,3), state.feasibility_region(2,3); ...
                                     state.feasibility_region(4,3), state.feasibility_region(3,3), state.feasibility_region(3,3), state.feasibility_region(4,3), state.feasibility_region(4,3)];
            
            feas_region_patch = patch(rectangle_feasibility(1,:), rectangle_feasibility(2,:), obj.color);
            feas_region_patch.FaceAlpha = 0.3;
            feas_region_patch.EdgeColor = [0.2660 0.3740 0.8880];     
            
            if strcmp(obj.sim_parameters.sim_type, 'obstacle')
                
                for i = 1 : obj.sim_parameters.obstacle_number
                    
                    x_m = obj.sim_parameters.obstacles(i, 1) - obj.sim_parameters.obstacles(i, 3) / 2;
                    x_M = obj.sim_parameters.obstacles(i, 1) + obj.sim_parameters.obstacles(i, 3) / 2;
                    y_m = obj.sim_parameters.obstacles(i, 2) - obj.sim_parameters.obstacles(i, 4) / 2; 
                    y_M = obj.sim_parameters.obstacles(i, 2) + obj.sim_parameters.obstacles(i, 4) / 2; 
                    object = [x_M, x_M, x_m, x_m, x_M; y_M, y_m, y_m, y_M, y_M];
                    object_patch = patch(object(1,:), object(2,:), 'k'); 
                    object_ptch.FaceAlpha = 0.3;
                    
                end
                
            end
            
%             feas_region_patch_estimate = patch(rectangle_feasibility_estimate(1,:), rectangle_feasibility_estimate(2,:), obj.color_estimate);
%             feas_region_patch_estimate.FaceAlpha = 0.1;
%             feas_region_patch_estimate.EdgeColor = [0.2660 0.3740 0.8880];            
            dcm = scatter(state.x(1,1) + state.x(2,1)/ obj.input.scheme_parameters.eta, ...
                        state.y(1,1) + state.y(2,1)/ obj.input.scheme_parameters.eta,...
                        'g','Linewidth',2);
            axis equal
            axis([-0.1 1.5 -0.25 0.55])   
            pbaspect([2 1 1]);
            legend([com, zmp, feas_region_patch, dcm], {'CoM', 'ZMP', 'current feas region', 'current dcm'});
            
        end
        
    end
    
    properties (Access = private)
    
        logs;
        sim_parameters;
        input;
        figure_handle;
        rectangle;
        initial_ds;
        color;
        color_estimate;
        counter;
        
    end
    
    
end
    
    