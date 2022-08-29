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
            
        end
        
        function [u, ftstp] = solve(obj, state, input)
        
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
                           + obj.tail_x;
                   x_u_M = (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) * obj.Delta_lambda / (obj.eta * obj.T_ds_0) ...
                           + obj.kinematic_and_temporal_buffer_x_M * obj.Delta_lambda + (obj.x_f_0 - obj.x_f_minus_one) * log(obj.Delta_lambda) / (obj.eta * obj.T_ds_0) ...        
                           + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                           + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_M) * exp( - obj.eta * obj.T_p) ...
                           - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                           + obj.tail_x;  
                   y_u_m = (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) * obj.Delta_lambda / (obj.eta * obj.T_ds_0) ...
                           + obj.kinematic_and_temporal_buffer_y_m * obj.Delta_lambda + (obj.y_f_0 - obj.y_f_minus_one) * log(obj.Delta_lambda) / (obj.eta * obj.T_ds_0) ...        
                           + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                           - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_m) * exp( - obj.eta * obj.T_p) ...
                           - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                           + obj.tail_y;
                   y_u_M = (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) * obj.Delta_lambda / (obj.eta * obj.T_ds_0) ...
                           + obj.kinematic_and_temporal_buffer_y_M * obj.Delta_lambda + (obj.y_f_0 - obj.y_f_minus_one) * log(obj.Delta_lambda) / (obj.eta * obj.T_ds_0) ...        
                           + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                           + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_M) * exp( - obj.eta * obj.T_p) ...
                           - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                           + obj.tail_y; 
                       
                   obj.mean = exp( - obj.eta * obj.T_ss_0) - (exp( - obj.eta * obj.T_ss_0) - exp( - obj.eta * obj.T_s_0)) / 2;
                   
                   if (obj.x_f_0 - obj.x_f_minus_one) >= 0
                                           
                       obj.linear_feasibility_region.a_x_m = obj.kinematic_and_temporal_buffer_x_m + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                              (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                       obj.linear_feasibility_region.b_x_m = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_m) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                              - ((obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0));
                                                          
                       obj.linear_feasibility_region.a_x_M = obj.kinematic_and_temporal_buffer_x_M + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                             - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                             + (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                       obj.linear_feasibility_region.b_x_M = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_M) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                              - ((obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                              / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0);
                                                                                                                          
                   else

                       obj.linear_feasibility_region.a_x_m = obj.kinematic_and_temporal_buffer_x_m + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                             - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                             + (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                       obj.linear_feasibility_region.b_x_m = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_m) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                              - ((obj.x_f_0 - obj.x_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                              / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0);                       
                       
                       obj.linear_feasibility_region.a_x_M = obj.kinematic_and_temporal_buffer_x_M + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                              (obj.x_f_minus_one - obj.x_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                       obj.linear_feasibility_region.b_x_M = + (obj.x_f_0 - obj.x_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.x_f_minus_one + (obj.kinematic_buffer_x_M) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.x_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_x + ( (obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                              - ((obj.x_f_0 - obj.x_f_minus_one) / ( obj.eta * obj.T_ds_0));
                                                                                                                                                                                 
                   end
                   
                   if (obj.y_f_0 - obj.y_f_minus_one) >= 0
                       
                       obj.linear_feasibility_region.a_y_m = obj.kinematic_and_temporal_buffer_y_m + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                              (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                       obj.linear_feasibility_region.b_y_m = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_m) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                              - ((obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0));
                                                          
                       obj.linear_feasibility_region.a_y_M = obj.kinematic_and_temporal_buffer_y_M + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                             - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                             + (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                       obj.linear_feasibility_region.b_y_M = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_M) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                              - ((obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                              / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0);
                       
                   else
                       
                       obj.linear_feasibility_region.a_y_m = obj.kinematic_and_temporal_buffer_y_m + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * (log( exp( - obj.eta * obj.T_ss_0)) ...
                                                             - log(exp(- obj.eta * obj.T_s_0))) / (exp( - obj.eta * obj.T_ss_0 ) - exp( - obj.eta * obj.T_s_0 )) ...
                                                             + (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);                                                                                                                 
                       obj.linear_feasibility_region.b_y_m = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              - obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_m) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0) )*(log( exp( - obj.eta * obj.T_s_0 ))) ...
                                                              - ((obj.y_f_0 - obj.y_f_minus_one) / (obj.eta * obj.T_ds_0)) * (log(exp(-obj.eta * obj.T_ss_0)) - log(exp(-obj.eta * obj.T_s_0))) ...
                                                              / (exp(-obj.eta * obj.T_ss_0) - exp(-obj.eta * obj.T_s_0)) * exp(- obj.eta * obj.T_s_0);                       
                       
                       obj.linear_feasibility_region.a_y_M = obj.kinematic_and_temporal_buffer_y_M + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) ) * ( 1 / obj.mean ) + ...
                                                              (obj.y_f_minus_one - obj.y_f_0) * exp(obj.eta * obj.T_ss_0) / (obj.eta * obj.T_ds_0);
                       obj.linear_feasibility_region.b_y_M = + (obj.y_f_0 - obj.y_f_minus_one) * (1 + obj.eta * obj.T_s_0) / (obj.eta * obj.T_ds_0) ...
                                                              + obj.d_z * (1 - exp(-obj.eta * obj.T_c)) / 2 + obj.y_f_minus_one + (obj.kinematic_buffer_y_M) * exp( - obj.eta * obj.T_p) ...
                                                              - obj.y_f_0 * exp( - obj.eta * obj.T_c) ...
                                                              + obj.tail_y + ( (obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0) )*( log(obj.mean) ) ...
                                                              - ((obj.y_f_0 - obj.y_f_minus_one) / ( obj.eta * obj.T_ds_0));
                       
                   end
                   
               else
                   
                   x_u_m = + obj.kinematic_and_temporal_buffer_x_m * obj.Delta_lambda ...
                           + (obj.x_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_x_m * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_x;
                   x_u_M = + obj.kinematic_and_temporal_buffer_x_M * obj.Delta_lambda ...
                           + (obj.x_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_x_M * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_x;
                   y_u_m = + obj.kinematic_and_temporal_buffer_y_m * obj.Delta_lambda ...
                           + (obj.y_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_y_m * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_y;
                   y_u_M = + obj.kinematic_and_temporal_buffer_y_M * obj.Delta_lambda ...
                           + (obj.y_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) ...
                           + obj.kinematic_buffer_y_M * exp( - obj.eta * obj.T_p) ...
                           + obj.tail_y;   
                       
                   obj.linear_feasibility_region.a_x_m = obj.kinematic_and_temporal_buffer_x_m;
                   obj.linear_feasibility_region.b_x_m = (obj.x_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_x_m * exp( - obj.eta * obj.T_p) + obj.tail_x;
                   obj.linear_feasibility_region.a_x_M = obj.kinematic_and_temporal_buffer_x_M;
                   obj.linear_feasibility_region.b_x_M = (obj.x_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_x_M * exp( - obj.eta * obj.T_p) + obj.tail_x;
                   obj.linear_feasibility_region.a_y_m = obj.kinematic_and_temporal_buffer_y_m;
                   obj.linear_feasibility_region.b_y_m = (obj.y_f_minus_one - obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_y_m * exp( - obj.eta * obj.T_p) + obj.tail_y;
                   obj.linear_feasibility_region.a_y_M = obj.kinematic_and_temporal_buffer_y_M;
                   obj.linear_feasibility_region.b_y_M = (obj.y_f_minus_one + obj.d_z / 2) * (1 - exp( - obj.eta * obj.T_c)) + obj.kinematic_buffer_y_M * exp( - obj.eta * obj.T_p) + obj.tail_y;

               end
               
               x_u_m_bar = obj.linear_feasibility_region.a_x_m * obj.Delta_lambda + obj.linear_feasibility_region.b_x_m;
               x_u_M_bar = obj.linear_feasibility_region.a_x_M * obj.Delta_lambda + obj.linear_feasibility_region.b_x_M;
               y_u_m_bar = obj.linear_feasibility_region.a_y_m * obj.Delta_lambda + obj.linear_feasibility_region.b_y_m;
               y_u_M_bar = obj.linear_feasibility_region.a_y_M * obj.Delta_lambda + obj.linear_feasibility_region.b_y_M;
               obj.feasibility_region(:, i) = [x_u_m; x_u_M; y_u_m; y_u_M; x_u_m_bar; x_u_M_bar; y_u_m_bar; y_u_M_bar];   
               
           end
           
        end

        function result = adaptTiming(obj)
            result = 0;
        end

        function result = selectKinematicConstraint(obj)
            result = 0;
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
        
    end
            
end