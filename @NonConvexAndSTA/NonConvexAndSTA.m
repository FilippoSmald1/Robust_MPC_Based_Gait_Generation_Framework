classdef NonConvexAndSTA < FeasibilityDrivenBase & handle
    
    methods (Access = public) 
        
        function obj = NonConvexAndSTA(input)
            
            obj.input = input;
            obj.feasibility_region = zeros(8, input.kar.number_of_subregions);
            obj.feasibility_region_linear_estimate = zeros(4, 1);
            obj.centerline_multiplier = obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * ...
                                        exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (0 : obj.input.scheme_parameters.C - 1));
            obj.tail_multiplier = obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * ...
                                        exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta * (obj.input.scheme_parameters.C : obj.input.scheme_parameters.P - 1));            
            obj.eta = obj.input.scheme_parameters.eta;
            obj.T_c = obj.input.scheme_parameters.T_c;
            obj.T_p = obj.input.scheme_parameters.T_p;
            obj.d_z = obj.input.scheme_parameters.d_z;
            obj.M = obj.input.scheme_parameters.M;   
            obj.T_ds_0 = input.footstep_plan.ds_duration;
            obj.lambda_j = zeros(obj.M, 1);
            obj.ni_j = zeros(obj.M, 1);
            obj.linear_feasibility_region = struct;
            obj.linear_feasibility_region.a_x_m = 0;
            obj.linear_feasibility_region.b_x_m = 0;
            obj.linear_feasibility_region.a_x_M = 0;
            obj.linear_feasibility_region.b_x_M = 0;
            obj.linear_feasibility_region.a_y_m = 0;
            obj.linear_feasibility_region.b_y_m = 0;
            obj.linear_feasibility_region.a_y_M = 0;
            obj.linear_feasibility_region.b_y_M = 0;
            obj.H_timing_qp = blkdiag(1, 10000000, 10000000);
            obj.f_timing_qp = zeros(3, 1);
            obj.A_timing_qp = zeros(7, 3);
            obj.b_timing_qp = zeros(7, 1);
            obj.t_min = 0.2;
            obj.t_max = 1.5;
            obj.eps_margin = 0.01;
            obj.eps_t = 0;
            obj.total_time_modification = 0;
                    
        end
        
        function [result, kin_constr] = solve(obj, state, input)

            obj.index = state.world_time_iter - round( obj.input.footstep_plan.timings(state.footstep_counter, 1) / obj.input.scheme_parameters.delta ) + 1;
            obj.input = input;
            
            if obj.index <= 1
                obj.total_time_modification
                obj.total_time_modification = 0;
            end
            computeFeasibilityRegion(obj, state, input);
            kin_constr = selectKinematicConstraint(obj, state);
            result = adaptTiming(obj, state);
            
        end
        
        function result = getFeasibilityRegion(obj)
            
            result = obj.feasibility_region;
            
        end
        
        function obj = computeFeasibilityRegion(obj, state, input)
           
           obj.input = input;
           obj.index = state.world_time_iter - round( obj.input.footstep_plan.timings(state.footstep_counter, 1) / obj.input.scheme_parameters.delta ) + 1;           
           obj.Delta_lambda = exp( - obj.eta * (obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) - state.world_time_iter * obj.input.scheme_parameters.delta));
           obj.T_s_0 = obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) - obj.input.footstep_plan.timings(state.footstep_counter, 1);
           obj.T_ss_0 = obj.T_s_0 - obj.T_ds_0;
           obj.tail_x = obj.tail_multiplier * obj.input.footstep_plan.tail_x + obj.input.footstep_plan.tail_x(end,1) * exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_p);
           obj.tail_y = obj.tail_multiplier * obj.input.footstep_plan.tail_y + obj.input.footstep_plan.tail_y(end,1) * exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_p);     
           
           for j = 1 : obj.M
               
              obj.lambda_j(j, 1) = exp( - obj.eta * (obj.input.footstep_plan.timings(state.footstep_counter + 1 + j, 1) - obj.input.footstep_plan.timings(state.footstep_counter + j, 1)));
              obj.ni_j(j, 1) = obj.T_ds_0 / (obj.input.footstep_plan.timings(state.footstep_counter + 1 + j, 1) - obj.input.footstep_plan.timings(state.footstep_counter + j, 1));
                          
           end
            
           
           for i = 1 : input.kar.number_of_subregions
               
               x_u_m = 0;
               x_u_M = 0;
               y_u_m = 0;
               y_u_M = 0;
               x_u_m_bar = 0;
               x_u_M_bar = 0;
               y_u_m_bar = 0;
               y_u_M_bar = 0;               
               
               obj.kinematic_buffer_x_m = 0;
               obj.kinematic_and_temporal_buffer_x_m = 0;
               obj.kinematic_buffer_x_M = 0;
               obj.kinematic_and_temporal_buffer_x_M = 0;
               obj.kinematic_buffer_y_m = 0;
               obj.kinematic_and_temporal_buffer_y_m = 0;
               obj.kinematic_buffer_y_M = 0;
               obj.kinematic_and_temporal_buffer_y_M = 0;
               
               % get values
               obj.x_f_0 = obj.input.footstep_plan.positions(state.footstep_counter_rm, 1);
               obj.x_f_minus_one = obj.input.footstep_plan.positions(state.footstep_counter_rm - 1, 1);
               obj.y_f_0 = obj.input.footstep_plan.positions(state.footstep_counter_rm, 2);
               obj.y_f_minus_one = obj.input.footstep_plan.positions(state.footstep_counter_rm - 1, 2);
               if state.footstep_counter > 1
                   1;
               end
               
               if state.current_sf == "right" && obj.index <= obj.input.footstep_plan.ds_samples 
                   obj.sf_sign = - 1;
               end
               if state.current_sf == "right" && obj.index > obj.input.footstep_plan.ds_samples 
                   obj.sf_sign = + 1;
               end               
               if state.current_sf == "left" && obj.index <= obj.input.footstep_plan.ds_samples 
                   obj.sf_sign = + 1;
               end
               if state.current_sf == "left" && obj.index > obj.input.footstep_plan.ds_samples 
                   obj.sf_sign = - 1;
               end  
               
               for j = 1 : obj.M

                   if j == 1
                       
                       obj.kinematic_buffer_x_m = obj.kinematic_buffer_x_m + obj.input.kar.subregion_parameters(i, 1) / 2 - obj.input.kar.subregion_parameters(i, 2);
                       obj.kinematic_buffer_x_M = obj.kinematic_buffer_x_M - obj.input.kar.subregion_parameters(i, 1) / 2 - obj.input.kar.subregion_parameters(i, 2);
                       
                       obj.kinematic_and_temporal_buffer_x_m = obj.kinematic_and_temporal_buffer_x_m + ( + obj.input.kar.subregion_parameters(i, 1) / 2 - obj.input.kar.subregion_parameters(i, 2)) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));
                       obj.kinematic_and_temporal_buffer_x_M = obj.kinematic_and_temporal_buffer_x_M + ( - obj.input.kar.subregion_parameters(i, 1) / 2 - obj.input.kar.subregion_parameters(i, 2)) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));
                       
                       % alternate pm sign for y

                       obj.kinematic_buffer_y_m = obj.kinematic_buffer_y_m + obj.input.kar.subregion_parameters(i, 3) / 2 - obj.sf_sign * obj.input.kar.subregion_parameters(i, 4);
                       obj.kinematic_buffer_y_M = obj.kinematic_buffer_y_M - obj.input.kar.subregion_parameters(i, 3) / 2 - obj.sf_sign * obj.input.kar.subregion_parameters(i, 4);
                       
                       obj.kinematic_and_temporal_buffer_y_m = obj.kinematic_and_temporal_buffer_y_m + ( + obj.input.kar.subregion_parameters(i, 3) / 2 - obj.sf_sign * obj.input.kar.subregion_parameters(i, 4)) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));
                       obj.kinematic_and_temporal_buffer_y_M = obj.kinematic_and_temporal_buffer_y_M + ( - obj.input.kar.subregion_parameters(i, 3) / 2 - obj.sf_sign * obj.input.kar.subregion_parameters(i, 4)) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));
                       obj.sf_sign = - obj.sf_sign;
                                                           
                   else
                       
                       obj.kinematic_buffer_x_m = obj.kinematic_buffer_x_m + obj.input.scheme_parameters.d_ax_subsequent / 2;
                       obj.kinematic_buffer_x_M = obj.kinematic_buffer_x_M - obj.input.scheme_parameters.d_ax_subsequent / 2;
                       
                       obj.kinematic_and_temporal_buffer_x_m = obj.kinematic_and_temporal_buffer_x_m + ( + obj.input.kar.subregion_parameters(i, 1) / 2) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) * prod(obj.lambda_j(1 : j - 1, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));
                       obj.kinematic_and_temporal_buffer_x_M = obj.kinematic_and_temporal_buffer_x_M + ( - obj.input.kar.subregion_parameters(i, 1) / 2) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) * prod(obj.lambda_j(1 : j - 1, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));                       
                       
                       obj.kinematic_buffer_y_m = obj.kinematic_buffer_y_m + obj.input.scheme_parameters.d_ay_subsequent / 2 - obj.sf_sign * obj.input.scheme_parameters.ell_y_subsequent;
                       obj.kinematic_buffer_y_M = obj.kinematic_buffer_y_M - obj.input.scheme_parameters.d_ay_subsequent / 2 - obj.sf_sign * obj.input.scheme_parameters.ell_y_subsequent;
                       
                       obj.kinematic_and_temporal_buffer_y_m = obj.kinematic_and_temporal_buffer_y_m + ( + obj.input.scheme_parameters.d_ay_subsequent / 2 - obj.sf_sign * obj.input.scheme_parameters.ell_y_subsequent) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) * prod(obj.lambda_j(1 : j - 1, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));
                       obj.kinematic_and_temporal_buffer_y_M = obj.kinematic_and_temporal_buffer_y_M + ( - obj.input.scheme_parameters.d_ay_subsequent / 2 - obj.sf_sign * obj.input.scheme_parameters.ell_y_subsequent) * ...
                                                               (1 - obj.lambda_j(j, 1) ^ obj.ni_j(j, 1)) * prod(obj.lambda_j(1 : j - 1, 1)) / (obj.ni_j(j, 1) * log(obj.lambda_j(j, 1)));                                                         
                       obj.sf_sign = - obj.sf_sign;
                       
                   end
                   
               end
           
               if obj.index <= obj.input.footstep_plan.ds_samples
                        
                   x_u_m = (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) * obj.Delta_lambda / (obj.eta * obj.T_ds_0) ...
                           + obj.kinematic_and_temporal_buffer_x_m * obj.Delta_lambda + (obj.x_f_0 - obj.x_f_minus_one) * log(obj.Delta_lambda) / (obj.eta * obj.T_ds_0) ...        
                           + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                           - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_m) * exp( - obj.eta * obj.T_p) ...
                           - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                           + obj.tail_x - state.w_bar(1,1) / obj.eta ^ 2;
                   x_u_M = (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) * obj.Delta_lambda / (obj.eta * obj.T_ds_0) ...
                           + obj.kinematic_and_temporal_buffer_x_M * obj.Delta_lambda + (obj.x_f_0 - obj.x_f_minus_one) * log(obj.Delta_lambda) / (obj.eta * obj.T_ds_0) ...        
                           + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                           + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_M) * exp( - obj.eta * obj.T_p) ...
                           - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                           + obj.tail_x - state.w_bar(1,1) / obj.eta ^ 2;  
                   y_u_m = (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) * obj.Delta_lambda / (obj.eta * obj.T_ds_0) ...
                           + obj.kinematic_and_temporal_buffer_y_m * obj.Delta_lambda + (obj.y_f_0 - obj.y_f_minus_one) * log(obj.Delta_lambda) / (obj.eta * obj.T_ds_0) ...        
                           + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                           - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_m) * exp( - obj.eta * obj.T_p) ...
                           - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                           + obj.tail_y - state.w_bar(2,1) / obj.eta ^ 2;
                   y_u_M = (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) * obj.Delta_lambda / (obj.eta * obj.T_ds_0) ...
                           + obj.kinematic_and_temporal_buffer_y_M * obj.Delta_lambda + (obj.y_f_0 - obj.y_f_minus_one) * log(obj.Delta_lambda) / (obj.eta * obj.T_ds_0) ...        
                           + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                           + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_M) * exp( - obj.eta * obj.T_p) ...
                           - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                           + obj.tail_y - state.w_bar(2,1) / obj.eta ^ 2; 
                       
                   obj.mean = exp( - obj.eta * obj.T_ss_0) - (exp( - obj.eta * obj.T_ss_0) - exp( - obj.eta * obj.T_s_0)) / 2;
                   
                   if i == 1
                       
                       if (obj.x_f_0 - obj.x_f_minus_one) >= 0

                           obj.linear_feasibility_region.a_x_m = obj.kinematic_and_temporal_buffer_x_m + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                                  (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                           obj.linear_feasibility_region.b_x_m = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_m) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                                  - ((obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0)) - state.w_bar(1,1) / obj.eta ^ 2;

                           obj.linear_feasibility_region.a_x_M = obj.kinematic_and_temporal_buffer_x_M + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                                 - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                                 + (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                           obj.linear_feasibility_region.b_x_M = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_M) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                                  - ((obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                                  / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0) - state.w_bar(1,1) / obj.eta ^ 2;

                       else

                           obj.linear_feasibility_region.a_x_m = obj.kinematic_and_temporal_buffer_x_m + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                                 - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                                 + (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                           obj.linear_feasibility_region.b_x_m = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_m) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                                  - ((obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                                  / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0) - state.w_bar(1,1) / obj.eta ^ 2;                       

                           obj.linear_feasibility_region.a_x_M = obj.kinematic_and_temporal_buffer_x_M + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                                  (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                           obj.linear_feasibility_region.b_x_M = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_M) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                                  - ((obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0)) - state.w_bar(1,1) / obj.eta ^ 2;

                       end

                       if (obj.y_f_0 - obj.y_f_minus_one) >= 0

                           obj.linear_feasibility_region.a_y_m = obj.kinematic_and_temporal_buffer_y_m + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                                  (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                           obj.linear_feasibility_region.b_y_m = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_m) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                                  - ((obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0)) - state.w_bar(2,1) / obj.eta ^ 2;

                           obj.linear_feasibility_region.a_y_M = obj.kinematic_and_temporal_buffer_y_M + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                                 - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                                 + (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                           obj.linear_feasibility_region.b_y_M = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_M) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                                  - ((obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                                  / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0) - state.w_bar(2,1) / obj.eta ^ 2;

                       else

                           obj.linear_feasibility_region.a_y_m = obj.kinematic_and_temporal_buffer_y_m + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                                 - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                                 + (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                           obj.linear_feasibility_region.b_y_m = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_m) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                                  - ((obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                                  / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0) - state.w_bar(2,1) / obj.eta ^ 2;                       

                           obj.linear_feasibility_region.a_y_M = obj.kinematic_and_temporal_buffer_y_M + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                                  (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                           obj.linear_feasibility_region.b_y_M = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                                  + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_M) * exp( - obj.eta * obj.T_p) ...
                                                                  - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                                  + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                                  - ((obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0)) - state.w_bar(2,1) / obj.eta ^ 2;

                       end
                   
                   end
                   
                   
               else
                   
                   x_u_m = + obj.kinematic_and_temporal_buffer_x_m * obj.Delta_lambda ...
                           + (obj.x_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_x_m * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_x - state.w_bar(1,1) / obj.eta ^ 2;
                   x_u_M = + obj.kinematic_and_temporal_buffer_x_M * obj.Delta_lambda ...
                           + (obj.x_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_x_M * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_x - state.w_bar(1,1) / obj.eta ^ 2;
                   y_u_m = + obj.kinematic_and_temporal_buffer_y_m * obj.Delta_lambda ...
                           + (obj.y_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_y_m * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_y - state.w_bar(2,1) / obj.eta ^ 2;
                   y_u_M = + obj.kinematic_and_temporal_buffer_y_M * obj.Delta_lambda ...
                           + (obj.y_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_y_M * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_y - state.w_bar(2,1) / obj.eta ^ 2;   
                   
                   if i == 1
                       obj.linear_feasibility_region.a_x_m = obj.kinematic_and_temporal_buffer_x_m;
                       obj.linear_feasibility_region.b_x_m = (obj.x_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_x_m * exp( - obj.eta * obj.T_p) + obj.tail_x - state.w_bar(1,1) / obj.eta ^ 2;
                       obj.linear_feasibility_region.a_x_M = obj.kinematic_and_temporal_buffer_x_M;
                       obj.linear_feasibility_region.b_x_M = (obj.x_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_x_M * exp( - obj.eta * obj.T_p) + obj.tail_x - state.w_bar(1,1) / obj.eta ^ 2;
                       obj.linear_feasibility_region.a_y_m = obj.kinematic_and_temporal_buffer_y_m;
                       obj.linear_feasibility_region.b_y_m = (obj.y_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_y_m * exp( - obj.eta * obj.T_p) + obj.tail_y - state.w_bar(2,1) / obj.eta ^ 2;
                       obj.linear_feasibility_region.a_y_M = obj.kinematic_and_temporal_buffer_y_M;
                       obj.linear_feasibility_region.b_y_M = (obj.y_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_y_M * exp( - obj.eta * obj.T_p) + obj.tail_y - state.w_bar(2,1) / obj.eta ^ 2;
                   end
                   
               end
               
               x_u_m_bar = obj.linear_feasibility_region.a_x_m * obj.Delta_lambda + obj.linear_feasibility_region.b_x_m;
               x_u_M_bar = obj.linear_feasibility_region.a_x_M * obj.Delta_lambda + obj.linear_feasibility_region.b_x_M;
               y_u_m_bar = obj.linear_feasibility_region.a_y_m * obj.Delta_lambda + obj.linear_feasibility_region.b_y_m;
               y_u_M_bar = obj.linear_feasibility_region.a_y_M * obj.Delta_lambda + obj.linear_feasibility_region.b_y_M;
               obj.feasibility_region(:, i) = [x_u_m; x_u_M; y_u_m; y_u_M; x_u_m_bar; x_u_M_bar; y_u_m_bar; y_u_M_bar];   
               
           end
           
        end

        function result = adaptTiming(obj, state)
            
            obj.f_timing_qp(1, 1) = - obj.Delta_lambda;
            
            obj.A_timing_qp(1, :) = [obj.linear_feasibility_region.a_x_m, -1, 0];
            obj.b_timing_qp(1, 1) = state.x(1, 1) + state.x(2,1) / obj.eta - obj.linear_feasibility_region.b_x_m - obj.eps_margin;
            obj.A_timing_qp(2, :) = [-obj.linear_feasibility_region.a_x_M, 1, 0];
            obj.b_timing_qp(2, 1) = - state.x(1, 1) - state.x(2,1) / obj.eta + obj.linear_feasibility_region.b_x_M - obj.eps_margin;
            obj.A_timing_qp(3, :) = [obj.linear_feasibility_region.a_y_m, 0, -1];
            obj.b_timing_qp(3, 1) = state.y(1, 1) + state.y(2,1) / obj.eta - obj.linear_feasibility_region.b_y_m - obj.eps_margin;
            obj.A_timing_qp(4, :) = [-obj.linear_feasibility_region.a_y_M, 0, 1];
            obj.b_timing_qp(4, 1) = - state.y(1, 1) - state.y(2,1) / obj.eta + obj.linear_feasibility_region.b_y_M - obj.eps_margin;            
            
            obj.A_timing_qp(5, :) = [1, 0, 0];
            obj.b_timing_qp(5, 1) = exp( - obj.eta * (obj.t_min - obj.input.scheme_parameters.delta * obj.index));
            obj.A_timing_qp(6, :) = [-1, 0, 0];
            obj.b_timing_qp(6, 1) = -exp(- obj.eta * (obj.t_max - obj.input.scheme_parameters.delta * obj.index)); 
            
            if obj.index <= obj.input.footstep_plan.ds_samples
                obj.eps_t = obj.T_ss_0 + max( abs(obj.input.footstep_plan.zmp_centerline_x(1,1) - obj.x_f_0), abs(obj.input.footstep_plan.zmp_centerline_y(1,1) - obj.y_f_0) ) / obj.input.scheme_parameters.v_max;
            else
                % using the centerline instead of the swing foot velocity
                % as in Matlab we do not have moving feet
                obj.eps_t = max( abs(obj.input.footstep_plan.zmp_centerline_x(1,1) - obj.x_f_0), abs(obj.input.footstep_plan.zmp_centerline_y(1,1) - obj.y_f_0) ) / obj.input.scheme_parameters.v_max;                
            end
                
            obj.A_timing_qp(7, :) = [1, 0, 0];   
            obj.b_timing_qp(7, 1) = exp( - obj.eta * obj.eps_t);
            
            sta = quadprog(obj.H_timing_qp, obj.f_timing_qp, obj.A_timing_qp, obj.b_timing_qp);
            
            result = state.world_time_iter;
            if ~isempty(sta)
               if abs(sta(1, 1) - obj.Delta_lambda) > 0.0001 
                   if obj.index <= obj.input.footstep_plan.ds_samples
                       if obj.input.footstep_plan.ds_samples - obj.index > 6
                           obj.total_time_modification = obj.total_time_modification + (result - round((obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) - log(1 / sta(1, 1)) / obj.eta) / obj.input.scheme_parameters.delta));
                           result = round((obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) - log(1 / sta(1, 1)) / obj.eta) / obj.input.scheme_parameters.delta);    
                       end
                   else
                       if round( obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) / obj.input.scheme_parameters.delta ) - state.world_time_iter > 6
                           obj.total_time_modification = obj.total_time_modification + (result - round((obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) - log(1 / sta(1, 1)) / obj.eta) / obj.input.scheme_parameters.delta));
                           result = round((obj.input.footstep_plan.timings(state.footstep_counter + 1, 1) - log(1 / sta(1, 1)) / obj.eta) / obj.input.scheme_parameters.delta);    
                       end                   
                   end
               end     
            end
                      
        end

        function result = selectKinematicConstraint(obj, state)
            
            x_u = state.x(1, 1) + state.x(2, 1) / obj.eta;
            y_u = state.y(1, 1) + state.y(2, 1) / obj.eta;
            i = 1;
            while ~(obj.feasibility_region(1, i) <= x_u && obj.feasibility_region(2, i) >= x_u && obj.feasibility_region(3, i) <= y_u && obj.feasibility_region(4, i) >= y_u)   
                i = i + 1; 
                if i > obj.input.kar.number_of_subregions
                    % if none of the regions is feasible, reinitialize the
                    % counter and pick the first region in the list, which
                    % has the highest priority
                    i = 1;
                    break;
                end
            end
            result = obj.input.kar.subregion_parameters(i, :);
        end
        
    end

    properties (Access = private)
      
        input;
        feasibility_region;
        feasibility_region_linear_estimate;
        index;
        eta;
        d_z;
        T_c;
        T_p;
        M;
        centerline_multiplier;
        tail_multiplier;
        tail_x;
        tail_y;
        x_f_0;
        x_f_minus_one;
        y_f_0;
        y_f_minus_one; 
        T_ds_0;
        T_ss_0;
        T_s_0;
        lambda_j;
        ni_j;
        Delta_lambda;
        kinematic_buffer_x_m;
        kinematic_and_temporal_buffer_x_m;
        kinematic_buffer_x_M;
        kinematic_and_temporal_buffer_x_M;        
        kinematic_buffer_y_m;
        kinematic_and_temporal_buffer_y_m;   
        kinematic_buffer_y_M;
        kinematic_and_temporal_buffer_y_M;  
        sf_sign;
        mean;
        linear_feasibility_region;  
        H_timing_qp;
        f_timing_qp;
        A_timing_qp;
        b_timing_qp;
        t_min;
        t_max;
        eps_margin;
        eps_t;
        total_time_modification;
        new_timings;
        
    end
            
end