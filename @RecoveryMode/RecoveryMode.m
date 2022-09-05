classdef RecoveryMode < FeasibilityDrivenBase & handle
    
    methods (Access = public)
        
        function obj = RecoveryMode(input)
            
            obj.input = input;
            obj.x_u_m = 0;
            obj.x_u_M = 0;
            obj.y_u_m = 0;
            obj.y_u_M = 0;
            obj.feasibility_region = [obj.x_u_m; obj.x_u_M; obj.y_u_m; obj.y_u_M];
            obj.centerline_multiplier = obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * ...
                                        exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (0 : obj.input.scheme_parameters.C - 1));
            obj.tail_multiplier = obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * ...
                                        exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (obj.input.scheme_parameters.C : obj.input.scheme_parameters.P - 1));            
            
            obj.P_matrix = obj.input.scheme_parameters.delta * tril( ones(obj.input.scheme_parameters.C, obj.input.scheme_parameters.C) );

            obj.A_stab_constr = zeros(1, obj.input.scheme_parameters.C + obj.input.scheme_parameters.M);
            for i = 1 : obj.input.scheme_parameters.C
                 obj.A_stab_constr(1, i) = (1 / obj.input.scheme_parameters.eta) * ...
                                  ( exp( - (i-1) * obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) - ...
                                    exp( - i * obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) );
            end
            obj.b_stab_constr = zeros(1);
            
            obj.footstep_weight = obj.input.scheme_parameters.footstep_weight_in_cost_function;
            obj.zmp_weight = obj.input.scheme_parameters.zmp_track_in_cost_function;
            
            obj.H = zeros(obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, obj.input.scheme_parameters.C + obj.input.scheme_parameters.M);
            obj.f = zeros(obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1);
            obj.H(1 : obj.input.scheme_parameters.C, 1 : obj.input.scheme_parameters.C) = eye(obj.input.scheme_parameters.C, obj.input.scheme_parameters.C);
            obj.H(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M) = ...
                 obj.footstep_weight * eye(obj.input.scheme_parameters.M, obj.input.scheme_parameters.M);
                         
            obj.d_ax_subsequent = obj.input.scheme_parameters.d_ax_subsequent; % this is never modified
            obj.d_ay_subsequent = obj.input.scheme_parameters.d_ay_subsequent;
            obj.ell_y_subsequent = obj.input.scheme_parameters.ell_subsequent;      
            
            obj.kin_constr_diff = eye(obj.input.scheme_parameters.M);
            for i = 2:obj.input.scheme_parameters.M
                obj.kin_constr_diff(i,i-1) = - 1;
            end 
            obj.A_kinematic_constr = [zeros(obj.input.scheme_parameters.M, obj.input.scheme_parameters.C), obj.kin_constr_diff; ...
                                      zeros(obj.input.scheme_parameters.M, obj.input.scheme_parameters.C), -obj.kin_constr_diff];
            obj.b_kinematic_constr = zeros(2 * obj.input.scheme_parameters.M, 1);
            obj.kin_constr_buffer = zeros(obj.input.scheme_parameters.M, 1);
            obj.kin_constr_buffer_y = zeros(obj.input.scheme_parameters.M, 2);
            obj.kin_multiplier = zeros(obj.input.scheme_parameters.M, 1);
            obj.kin_multiplier(1, 1) = 1;
            
            obj.A_zmp_constr = zeros(2 * obj.input.scheme_parameters.C, obj.input.scheme_parameters.C + obj.input.scheme_parameters.M);
            obj.b_zmp_constr = zeros(2 * obj.input.scheme_parameters.C, 1);
            obj.A_zmp_constr(1 : obj.input.scheme_parameters.C, 1 : obj.input.scheme_parameters.C) = obj.P_matrix;
            obj.A_zmp_constr(obj.input.scheme_parameters.C + 1 : 2 * obj.input.scheme_parameters.C, 1 : obj.input.scheme_parameters.C) = - obj.P_matrix;

            obj.A_eq_ss = obj.A_stab_constr;
            obj.b_eq_ss = obj.b_stab_constr;
            obj.A_eq_ds = [obj.A_stab_constr; zeros(1, obj.input.scheme_parameters.C), 1, zeros(1, obj.input.scheme_parameters.M - 1)]; 
            obj.b_eq_ds = [obj.b_stab_constr; 0];             

            obj.A_ineq = [obj.A_zmp_constr; obj.A_kinematic_constr];
            obj.b_ineq = [obj.b_zmp_constr; obj.b_kinematic_constr]; 
            
            obj.options =  optimoptions(@quadprog, 'Display','off');
                       
        end
        
        function [u, ftstp] = solve(obj, state, input)
            
            obj.input = input;  
            obj.feasibility_region = zeros(8, obj.input.kar.number_of_subregions);
            
            obj.d_ax = obj.input.scheme_parameters.d_ax; % this is the default kinematic admissible region (can be modified)
            obj.d_ay = obj.input.scheme_parameters.d_ay;
            obj.ell_y = obj.input.scheme_parameters.ell_y;
            obj.ell_x = obj.input.scheme_parameters.ell_x;
            
            obj.A_zmp_constr(1 : obj.input.scheme_parameters.C, obj.input.scheme_parameters.C + 1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M) = - obj.input.footstep_plan.mapping(:, 2 : obj.input.scheme_parameters.M + 1);
            obj.A_zmp_constr(obj.input.scheme_parameters.C + 1 : 2 * obj.input.scheme_parameters.C, obj.input.scheme_parameters.C + 1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M) =  obj.input.footstep_plan.mapping(:, 2 : obj.input.scheme_parameters.M + 1);
            
            obj.H(1 : obj.input.scheme_parameters.C, 1 : obj.input.scheme_parameters.C) = eye(obj.input.scheme_parameters.C) + obj.zmp_weight * obj.P_matrix' * obj.P_matrix;
            obj.H(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M) = ...
                 obj.footstep_weight * eye(obj.input.scheme_parameters.M, obj.input.scheme_parameters.M);
            
            obj.index = state.world_time_iter - round( obj.input.footstep_plan.timings(state.footstep_counter, 1) / obj.input.scheme_parameters.delta ) + 1;
            
            % x component           
            obj.b_stab_constr = state.x(1,1) + state.x(2,1) / obj.input.scheme_parameters.eta ...
                                 - state.x(3,1) * (1 - exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_c)) ...
                                 - obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                                 - obj.input.footstep_plan.tail_x(end,1) * exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_p) ...
                                 + state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2;
            
            obj.kin_constr_buffer(1, 1) = obj.d_ax / 2.0 + obj.ell_x;
            for i = 2 : obj.input.scheme_parameters.M
                obj.kin_constr_buffer(i, 1) = obj.d_ax_subsequent / 2.0;    
            end
            obj.b_kinematic_constr = [obj.kin_multiplier * state.sf_pos_ss(1,1) + obj.kin_constr_buffer; ...
                                      - obj.kin_multiplier * state.sf_pos_ss(1,1) + obj.kin_constr_buffer];
             
            obj.b_zmp_constr(1 : obj.input.scheme_parameters.C, 1) = - state.x(3,1) + obj.input.scheme_parameters.d_z / 2 + obj.input.footstep_plan.mapping(:, 1) * state.sf_pos_ss(1,1); 
            obj.b_zmp_constr(obj.input.scheme_parameters.C + 1 : 2 * obj.input.scheme_parameters.C, 1) = + state.x(3,1) + obj.input.scheme_parameters.d_z / 2 - obj.input.footstep_plan.mapping(:, 1) * state.sf_pos_ss(1,1);                                 
            
            obj.b_ineq = [obj.b_zmp_constr; obj.b_kinematic_constr];
            
            obj.f(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1) = ...
                    - obj.footstep_weight * obj.input.footstep_plan.positions(state.footstep_counter_rm : state.footstep_counter_rm + obj.input.scheme_parameters.M - 1, 1);           
            
            obj.f(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1) = ...
                    - obj.footstep_weight * obj.input.footstep_plan.positions(state.footstep_counter_rm : state.footstep_counter_rm + obj.input.scheme_parameters.M - 1, 1);
                                     
            obj.f(1 : obj.input.scheme_parameters.C, 1) = ...                   
                    -  obj.zmp_weight * obj.input.footstep_plan.zmp_centerline_x;
               
                
            if obj.index <= obj.input.footstep_plan.ds_samples 
               obj.A_ineq = obj.A_zmp_constr;
               obj.b_ineq = obj.b_zmp_constr; 
               obj.b_eq_ds(1,1) = obj.b_stab_constr;
               obj.b_eq_ds(2,1) = obj.input.footstep_plan.positions(state.footstep_counter_rm, 1);           
               solution = quadprog(obj.H, obj.f, obj.A_ineq, obj.b_ineq, obj.A_eq_ds, obj.b_eq_ds, [], [], [], obj.options);
            else
               obj.A_ineq = [obj.A_zmp_constr; obj.A_kinematic_constr];
               obj.b_ineq = [obj.b_zmp_constr; obj.b_kinematic_constr];                
               obj.A_eq_ss = obj.A_stab_constr;   
               obj.b_eq_ss = obj.b_stab_constr;
%                if round( obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) / obj.input.scheme_parameters.delta ) - state.world_time_iter <= 2 
%                   obj.A_eq_ss = [obj.A_eq_ss; zeros(1, obj.input.scheme_parameters.C + obj.input.scheme_parameters.M)];   
%                   obj.A_eq_ss(2, obj.input.scheme_parameters.C + 1) = 1;
%                   obj.b_eq_ss = [obj.b_eq_ss; obj.previous_ftstp(1, 1)];                   
%                end               
               solution = quadprog(obj.H, obj.f, obj.A_ineq, obj.b_ineq, obj.A_eq_ss, obj.b_eq_ss, [], [], [], obj.options);                
            end
            if isempty(solution)
                stop = 'stop';
            end
            u(1,1) = solution(1);
            ftstp(1,1) = solution(obj.input.scheme_parameters.C + 1);
            
             
            % y component
            obj.b_stab_constr = state.y(1,1) + state.y(2,1) / obj.input.scheme_parameters.eta ...
                                - state.y(3,1) * (1 - exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_c)) ...
                                - obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                                - obj.input.footstep_plan.tail_y(end,1) * exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_p) ...
                                + state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2;

            if state.current_sf == "right"
                obj.sf_sign = + 1;
            else
                obj.sf_sign = - 1;
            end

            obj.kin_constr_buffer_y(1, 1) = obj.d_ay / 2.0 + obj.sf_sign * obj.ell_y;
            obj.kin_constr_buffer_y(1, 2) = - obj.d_ay / 2.0 + obj.sf_sign * obj.ell_y;                

            for i = 2 : obj.input.scheme_parameters.M
                obj.sf_sign = - obj.sf_sign;
                obj.kin_constr_buffer_y(i, 1) = obj.d_ay_subsequent / 2.0 + obj.sf_sign * obj.ell_y_subsequent;  
                obj.kin_constr_buffer_y(i, 2) = - obj.d_ay_subsequent / 2.0 + obj.sf_sign * obj.ell_y_subsequent;
            end
            obj.b_kinematic_constr = [obj.kin_multiplier * state.sf_pos_ss(2,1) + obj.kin_constr_buffer_y(:, 1); ...
                                      - obj.kin_multiplier * state.sf_pos_ss(2,1) - obj.kin_constr_buffer_y(:, 2)];                            

            obj.b_zmp_constr(1 : obj.input.scheme_parameters.C, 1) = - state.y(3,1) + obj.input.scheme_parameters.d_z / 2 + obj.input.footstep_plan.mapping(:, 1) * state.sf_pos_ss(2,1); 
            obj.b_zmp_constr(obj.input.scheme_parameters.C + 1 : 2 * obj.input.scheme_parameters.C, 1) = + state.y(3,1) + obj.input.scheme_parameters.d_z / 2 - obj.input.footstep_plan.mapping(:, 1) * state.sf_pos_ss(2,1);                                                 
          
            obj.f(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1) = ...
                    - obj.footstep_weight * obj.input.footstep_plan.positions(state.footstep_counter_rm : state.footstep_counter_rm + obj.input.scheme_parameters.M - 1, 2);

            obj.f(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1) = ...
                    - obj.footstep_weight * obj.input.footstep_plan.positions(state.footstep_counter_rm : state.footstep_counter_rm + obj.input.scheme_parameters.M - 1, 2);
                                     
            obj.f(1 : obj.input.scheme_parameters.C, 1) = ...                   
                    -  obj.zmp_weight * obj.input.footstep_plan.zmp_centerline_y;                
                
             if obj.index <= obj.input.footstep_plan.ds_samples 
               obj.A_ineq = obj.A_zmp_constr;
               obj.b_ineq = obj.b_zmp_constr;
               obj.b_eq_ds(1,1) = obj.b_stab_constr;
               obj.b_eq_ds(2,1) = obj.input.footstep_plan.positions(state.footstep_counter_rm, 2);    
               solution = quadprog(obj.H, obj.f, obj.A_ineq, obj.b_ineq, obj.A_eq_ds, obj.b_eq_ds, [], [], [], obj.options);
             else
               obj.A_ineq = [obj.A_zmp_constr; obj.A_kinematic_constr];
               obj.b_ineq = [obj.b_zmp_constr; obj.b_kinematic_constr];                 
               obj.A_eq_ss = obj.A_stab_constr;   
               obj.b_eq_ss = obj.b_stab_constr;
%                if round( obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) / obj.input.scheme_parameters.delta ) - state.world_time_iter <= 1  
%                   obj.A_eq_ss = [obj.A_eq_ss; zeros(1, obj.input.scheme_parameters.C + obj.input.scheme_parameters.M)];   
%                   obj.A_eq_ss(2, obj.input.scheme_parameters.C + 1) = 1;
%                   obj.b_eq_ss = [obj.b_eq_ss; obj.previous_ftstp(2, 1)];                   
%                end
               solution = quadprog(obj.H, obj.f, obj.A_ineq, obj.b_ineq, obj.A_eq_ss, obj.b_eq_ss, [], [], [], obj.options);               
             end
            if isempty(solution)
                stop = 'stop';
            end            
            u(2,1) = solution(1);
            ftstp(2,1) = solution(obj.input.scheme_parameters.C + 1,1);
            ftstp(3,1) = 0;
            obj.previous_ftstp = ftstp; 
            
        end
        
        function result = getFeasibilityRegion(obj)
            result = 0;
            
        end
        
        function obj = computeFeasibilityRegion(obj, state, input)
            result = 0;
            % put also here the computation!!!
        end
        
    end
    
    properties (Access = private)
        
        input;
        x_u_m;
        x_u_M;
        y_u_m;
        y_u_M;
        x_u;
        y_u;
        feasibility_region;
        centerline_multiplier;
        tail_multiplier;
        P_matrix;
        H;
        f;
        A_eq_ss;
        b_eq_ss;
        A_eq_ds;
        b_eq_ds;
        A_stab_constr;
        b_stab_constr;
        A_ds_constr;
        b_ds_constr;
        A_ineq;
        b_ineq;
        A_zmp_constr;
        b_zmp_constr;
        A_kinematic_constr;
        b_kinematic_constr;
        d_ax;
        d_ay;
        ell_y;
        ell_x;
        d_ax_subsequent;
        d_ay_subsequent;
        ell_y_subsequent; 
        kin_constr_diff;
        kin_constr_buffer;
        kin_multiplier;
        kin_constr_buffer_y;
        sf_sign;
        footstep_weight;
        zmp_weight;
        index;
        previous_ftstp;
        options;
                   
    end
    
end