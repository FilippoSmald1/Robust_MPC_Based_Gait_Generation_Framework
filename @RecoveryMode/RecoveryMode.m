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

            obj.A_stab_constr = zeros(1, obj.input.scheme_parameters.C);
            for i = 1 : obj.input.scheme_parameters.C
                 obj.A_stab_constr(1, i) = (1 / obj.input.scheme_parameters.eta) * ...
                                  ( exp( - (i-1) * obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) - ...
                                    exp( - i * obj.input.scheme_parameters.eta * obj.input.scheme_parameters.delta) );
            end
            obj.b_stab_constr = zeros(1);
            obj.H = zeros(obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, obj.input.scheme_parameters.C + obj.input.scheme_parameters.M);
            obj.f = zeros(obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1);
            obj.H(1 : obj.input.scheme_parameters.C, 1 : obj.input.scheme_parameters.C) = eye(obj.input.scheme_parameters.C, obj.input.scheme_parameters.C);
            obj.H(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M) = ...
                 eye(obj.input.scheme_parameters.M, obj.input.scheme_parameters.M);
            obj.mapping = zeros(obj.input.scheme_parameters.C, obj.input.scheme_parameters.M + 1);
             
            obj.d_ax = obj.input.scheme_parameters.d_ax; % this is the default kinematic admissible region (can be modified)
            obj.d_ay = obj.input.scheme_parameters.d_ay;
            obj.ell = obj.input.scheme_parameters.ell;
            obj.d_ax_subsequent = obj.input.scheme_parameters.d_ax_subsequent; % this is never modified
            obj.d_ay_subsequent = obj.input.scheme_parameters.d_ay_subsequent;
            obj.ell_subsequent = obj.input.scheme_parameters.ell_subsequent;            
            
        end
        
        function [u, ftstp] = solve(obj, state, input)
            
            obj.input = input;    
            
            % x component           
             obj.b_stab_constr = state.x(1,1) + state.x(2,1) / obj.input.scheme_parameters.eta ...
                                 - state.x(3,1) * (1 - exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_c)) ...
                                 - obj.tail_multiplier * obj.input.footstep_plan.tail_x ...
                                 + state.w_bar(1,1) / obj.input.scheme_parameters.eta ^ 2;
                             
            obj.f(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1) = ...
                 - 2 * obj.input.footstep_counter(state.footstep_counter + 1 : state.footstep_counter + 1 + obj.input.scheme_parameters.M, 1);
                   
            % y component
            obj.b_stab_constr = state.y(1,1) + state.y(2,1) / obj.input.scheme_parameters.eta ...
                                - state.y(3,1) * (1 - exp( - obj.input.scheme_parameters.eta * obj.input.scheme_parameters.T_c)) ...
                                - obj.tail_multiplier * obj.input.footstep_plan.tail_y ...
                                + state.w_bar(2,1) / obj.input.scheme_parameters.eta ^ 2;
                   
            obj.f(obj.input.scheme_parameters.C+1 : obj.input.scheme_parameters.C + obj.input.scheme_parameters.M, 1) = ...
                 - 2 * obj.input.footstep_counter(state.footstep_counter + 1 : state.footstep_counter + 1 + obj.input.scheme_parameters.M, 2);

        end
        
        function result = getFeasibilityRegion(obj)
            result = 0;
        end
        
        function obj = computeFeasibilityRegion(obj, state, input)
            result = 0;
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
        A_eq;
        b_eq;
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
        ell;
        d_ax_subsequent;
        d_ay_subsequent;
        ell_subsequent;        
                   
    end
    
end