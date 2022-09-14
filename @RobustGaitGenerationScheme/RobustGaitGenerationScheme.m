classdef RobustGaitGenerationScheme < handle

    methods (Access = public)

        function obj = RobustGaitGenerationScheme(input, state, sim_parameters)
            
            % constructor
            obj.input = input;
            obj.sim_parameters = sim_parameters;
            obj.sm_instance = StandardMode(obj.input);
            obj.rm_instance = RecoveryMode(obj.input);
            obj.rm_feasibility_processing = NonConvexAndSTA(obj.input);
            obj.dob_instance = DisturbanceObserver(obj.input, state);
            
            obj.steps_in_horizon = zeros(obj.input.scheme_parameters.F + 2, 2);   
            obj.centerline_temp_temp_x = zeros(2 * obj.input.scheme_parameters.P, 1);
            obj.centerline_temp_temp_y = zeros(2 * obj.input.scheme_parameters.P, 1); 
            obj.centerline_temp_x = cell(obj.input.scheme_parameters.F + 1, 1);
            obj.centerline_temp_y = cell(obj.input.scheme_parameters.F + 1, 1);
            obj.mapping_buffer = zeros(2 * obj.input.scheme_parameters.P, obj.input.scheme_parameters.F + 2);
            
            obj.u = [0; 0]; % zmp velocity command (x,y)
            obj.ftstp = [0; 0; 0];  % footstep position (x,y)
            
            ch = cosh(obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta);
            sh = sinh(obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta);
            eta = obj.input.scheme_parameters.eta; % temp variable to enlight the writing
            delta = obj.input.scheme_parameters.delta; % temp variable to enlight the writing
            obj.A_upd = [ch, sh/eta, 1-ch; 
                         eta*sh, ch, -eta*sh; 
                         0, 0, 1];
            obj.B_upd = [delta-sh/eta; 
                         1-ch; 
                         delta];
            obj.D_upd = [0;
                         delta;
                         0];
                     
            obj.mode = 'standard_mode';
            obj.obstacles_are_present = false;
            obj.feasibility_region = [0; 0; 0; 0];
            obj.new_timings = obj.input.footstep_plan.timings;
            
        end

        function state_ = update(obj, state)
            
             % update current footstep counter when new step begins
             if state.world_time_iter == round( (obj.input.footstep_plan.timings(state.footstep_counter, 1))/ obj.input.scheme_parameters.delta) + obj.input.footstep_plan.ds_samples
                state.sf_pos = state.next_sf_pos;
                if state.current_sf == "right"
                   state.current_sf = "left"; 
                else
                   state.current_sf = "right"; 
                end                
                state.footstep_counter_rm = state.footstep_counter_rm + 1;
                state.sf_pos_ss = state.next_sf_pos;  
                obj.input.footstep_plan.positions(state.footstep_counter, :) = state.sf_pos';            
                obj.mode = 'standard_mode';               
             end
             
             if state.world_time_iter == round(obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) / obj.input.scheme_parameters.delta)            
                state.step_time_iter = 1;           
                state.sf_pos = state.next_sf_pos; 
                state.sf_pos_ss = state.next_sf_pos;
                state.footstep_counter = state.footstep_counter + 1;
                state.footstep_counter_sm = state.footstep_counter + 1;                
                obj.new_timings(state.footstep_counter, 1) = state.sim_iter;                
                difference = obj.input.footstep_plan.positions(state.footstep_counter + 1, :) - state.sf_pos_ss';                 
                obj.input.footstep_plan.positions(state.footstep_counter + 1, :) = state.sf_pos';
                obj.input.footstep_plan.positions(state.footstep_counter + 2 : end, 1) = obj.input.footstep_plan.positions(state.footstep_counter + 2 : end, 1) - difference(1);
                obj.input.footstep_plan.positions(state.footstep_counter + 2 : end, 2) = obj.input.footstep_plan.positions(state.footstep_counter + 2 : end, 2) - difference(2);
                obj.input.footstep_plan.positions(state.footstep_counter + 2 : end, 3) = obj.input.footstep_plan.positions(state.footstep_counter + 2 : end, 3) - difference(3);                                
                obj.mode = 'standard_mode';                
             end
             
             % observer
             obj.dob_instance.update([state.x(1,1); state.x(3,1)], ...
                                     [state.y(1,1); state.y(3,1)], ...
                                     obj.u);  
             state.w_bar = obj.dob_instance.getDisturbance();
             
             
             % build a ZMP trajectory to be used when needed                                 
             obj.zmpTrajectoryGenerator(state);
             
             % feasibility check (only if cleared by the obstacle
             % detection)
             if obj.sm_instance.feasibilityCheck(state, obj.input) && ~obj.obstacles_are_present && strcmp(obj.mode, 'standard_mode') 
                obj.mode = 'standard_mode'; 
                if state.footstep_counter == 3 && strcmp(obj.sim_parameters.sim_type, 'obstacle')
                    obj.mode = 'recovery_mode';    
                end
             else
                obj.mode = 'recovery_mode'; 
             end
                   
             % control cycle
             if strcmp(obj.mode, 'standard_mode')
                 [obj.u, obj.ftstp] = obj.sm_instance.solve(state, obj.input);
                 obj.feasibility_region = obj.sm_instance.getFeasibilityRegion();
                 state.feasibility_region = obj.feasibility_region;
             else
                 obj.input.kar.number_of_subregions = 3;
                 obj.input.kar.subregion_parameters = [0.5, 0, 0.25, 0.2; ...
                                                       0.15, 0.175, 0.15, 0.00; ...
                                                       0.15, -0.175, 0.15, 0.00];                 
                 % obstacle detection 
                 if state.footstep_counter == 3 && strcmp(obj.sim_parameters.sim_type, 'obstacle')
                    obj.input.kar.number_of_subregions = 3;
                    obj.input.kar.subregion_parameters = [0.44, -0.03, 0.15, 0.15; ...
                                                          0.15, 0.175, 0.15, 0.00; ...
                                                          0.15, -0.175, 0.15, 0.00];
                 end
                 
                 [output, kin_constr] = obj.rm_feasibility_processing.solve(state, obj.input);               
                 if state.footstep_counter > 1 
                     state.world_time_iter = output;
                 end                
                 obj.input.scheme_parameters.d_ax = kin_constr(1, 1);    
                 obj.input.scheme_parameters.ell_x = kin_constr(1, 2); 
                 obj.input.scheme_parameters.d_ay = kin_constr(1, 3);
                 obj.input.scheme_parameters.ell_y = kin_constr(1, 4);  
                 [obj.u, obj.ftstp] = obj.rm_instance.solve(state, obj.input);                   
                 state.feasibility_region = obj.rm_feasibility_processing.getFeasibilityRegion(); 
             end
             
             % integrate LIP and update target swing foot target
             state = obj.integrateModel(state, obj.u);                      
             state.next_sf_pos = obj.ftstp;
                          
             % update counters
             state.world_time_iter = state.world_time_iter + 1;
             state.step_time_iter = state.step_time_iter + 1;
             
             % return updated state structure
             state_ = state; 
             
        end
        
        function init_state_proposal = proposeFeasibleInitialState(obj, state)
            
            zmpdot = zeros(2,1);
            obj.dob_instance.update([state.x(1,1); state.x(3,1)], ...
                                     [state.y(1,1); state.y(3,1)], ...
                                     zmpdot);  
            state.w_bar = obj.dob_instance.getDisturbance();
                         
            obj.zmpTrajectoryGenerator(state); 
            if ~obj.sm_instance.feasibilityCheck(state, obj.input)
                % if initial state is unfeasible propose a feasible initial
                % state (CoM position with zero velocity)
                region = obj.sm_instance.getFeasibilityRegion(); 
                init_state_proposal = [ (region(1,1) + region(2,1)) / 2; 0; ...
                                        (region(3,1) + region(4,1)) / 2; 0];                
            else
                % if initial state is feasible keep current state
                init_state_proposal = [ state.x(1,1); state.x(2,1); ...
                                        state.y(1,1); state.y(2,1)];                                    
            end

        end
        
        function obj = zmpTrajectoryGenerator(obj, state)
            
             index = state.world_time_iter - round( obj.input.footstep_plan.timings(state.footstep_counter, 1) / obj.input.scheme_parameters.delta ) + 1;
             obj.mapping_buffer(:,:) = 0;
                              
             obj.steps_in_horizon(1 : obj.input.scheme_parameters.F + 2, :) = ...
                 obj.input.footstep_plan.positions(state.footstep_counter : state.footstep_counter + obj.input.scheme_parameters.F + 1, 1:2);
             
             time_counter = 1;
             for i = 1 : obj.input.scheme_parameters.F + 1
                 
                 obj.ss_samples = round( (obj.input.footstep_plan.timings(i + 1, 1) - obj.input.footstep_plan.timings(i, 1)) / obj.input.scheme_parameters.delta ) - ...
                                  obj.input.footstep_plan.ds_samples;
                 obj.centerline_temp_x{i} = [obj.steps_in_horizon(i, 1) + (obj.steps_in_horizon(i + 1, 1) - obj.steps_in_horizon(i, 1)) * ...
                                            (1 : obj.input.footstep_plan.ds_samples )' / obj.input.footstep_plan.ds_samples; ...
                                            obj.steps_in_horizon(i + 1, 1) * ones(obj.ss_samples ,1)];
                 obj.centerline_temp_y{i} = [obj.steps_in_horizon(i, 2) + (obj.steps_in_horizon(i + 1, 2) - obj.steps_in_horizon(i, 2)) * ...
                                            (1 : obj.input.footstep_plan.ds_samples )' / obj.input.footstep_plan.ds_samples; ...
                                            obj.steps_in_horizon(i + 1, 2) * ones(obj.ss_samples ,1)]; 
                 
                 obj.centerline_temp_temp_x(time_counter : time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples - 1, 1) = obj.centerline_temp_x{i};                    
                 obj.centerline_temp_temp_y(time_counter : time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples - 1, 1) = obj.centerline_temp_y{i};                  
                 
                 if index <= obj.input.footstep_plan.ds_samples
                     obj.mapping_buffer(time_counter : time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples - 1, i + 1) = [(0 : obj.input.footstep_plan.ds_samples - 1)' / obj.input.footstep_plan.ds_samples; ...
                                                                                                                                        ones(obj.ss_samples ,1)];
                     obj.mapping_buffer(time_counter : time_counter + obj.input.footstep_plan.ds_samples - 1, i) = flip((1 : obj.input.footstep_plan.ds_samples )' / obj.input.footstep_plan.ds_samples);
                 else
                     obj.mapping_buffer(time_counter : time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples - 1, i) = [ones(obj.ss_samples ,1); ...
                                                                                                                                     flip((1 : obj.input.footstep_plan.ds_samples )' / obj.input.footstep_plan.ds_samples)];
                     obj.mapping_buffer(time_counter + obj.ss_samples : time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples - 1, i + 1) = (0 : obj.input.footstep_plan.ds_samples - 1)' / obj.input.footstep_plan.ds_samples;
                 end
                 
                 time_counter = time_counter + obj.ss_samples + obj.input.footstep_plan.ds_samples;
                 
                 
             end
             
             index = state.world_time_iter - round( obj.input.footstep_plan.timings(state.footstep_counter, 1) / obj.input.scheme_parameters.delta ) + 1;
             obj.input.footstep_plan.zmp_centerline_x = obj.centerline_temp_temp_x(index : index + obj.input.scheme_parameters.C - 1, 1);
             obj.input.footstep_plan.zmp_centerline_y = obj.centerline_temp_temp_y(index : index + obj.input.scheme_parameters.C - 1, 1);
             obj.input.footstep_plan.tail_x = obj.centerline_temp_temp_x(index + obj.input.scheme_parameters.C: index + obj.input.scheme_parameters.P - 1, 1);                                                                           
             obj.input.footstep_plan.tail_y = obj.centerline_temp_temp_y(index + obj.input.scheme_parameters.C: index + obj.input.scheme_parameters.P - 1, 1);                                                                                                
             if index <= obj.input.footstep_plan.ds_samples
                 obj.input.footstep_plan.mapping = obj.mapping_buffer(index : index + obj.input.scheme_parameters.C - 1, 1 : obj.input.scheme_parameters.M + 1);
             else
                 obj.input.footstep_plan.mapping = obj.mapping_buffer(index - obj.input.footstep_plan.ds_samples : index + obj.input.scheme_parameters.C - 1 - obj.input.footstep_plan.ds_samples, 1 : obj.input.scheme_parameters.M + 1);                 
             end
        end
        
        % getters
        
        function w_bar = getDisturbance(obj)

            w_bar = obj.dob_instance.getDisturbance();

        end
        
        function data = getModifiedInputStructure(obj)

            data = obj.input;

        end     
        
        function data = getCurrentFeasibilityRegion(obj)

            data = obj.input;

        end     
        
        function timings = getNewTimings(obj)
            
            timings = obj.new_timings;
            
        end

    end
    
    methods (Access = private) 
       
        function state_ = integrateModel(obj, state, u)
            
            state_ = state;
            state_.x = obj.A_upd * state.x + obj.B_upd * u(1,1);
            state_.y = obj.A_upd * state.y + obj.B_upd * u(2,1);
            
        end
        
    end

    properties (Access = private)
        
        % building blocks
        input;
        sm_instance;
        rm_instance;
        rm_feasibility_processing;
        dob_instance;
        restriction_builder;
        sim_parameters;
        
        % variables
        u;
        ftstp;
        A_upd;
        B_upd;
        D_upd;
        mode;
        feasibility_region;
        obstacles_are_present;       
        
        % buffers
        steps_in_horizon;
        centerline_temp_x;
        centerline_temp_y;
        centerline_temp_temp_x;
        centerline_temp_temp_y;      
        ss_samples;
        mapping_buffer;
        new_timings;
        
    end
    
end
