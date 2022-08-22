classdef DisturbanceObserver < handle

    methods (Access = public)

        function obj = DisturbanceObserver(input, init_state)

            obj.input = input;  

            ch = cosh(input.scheme_parameters.eta*input.scheme_parameters.delta);
            sh = sinh(input.scheme_parameters.eta*input.scheme_parameters.delta);
            eta = input.scheme_parameters.eta; % temp variable for clarity
            delta = input.scheme_parameters.delta; 
            obj.A_upd_obs = [ch, sh/eta, 1-ch, 0; 
                             eta*sh, ch, -eta*sh, delta; 
                             0, 0, 1, 0; 
                             0, 0, 0, 1];
            obj.B_upd_obs = [delta-sh/eta; 
                             1-ch; 
                             delta; 
                             0];
            obj.C_measurable_output = [1 0 0 0; 
                                       0 0 1 0];
            p = [0.8, 0.7, 0.6, 0.5];
            obj.G = place(obj.A_upd_obs', obj.C_measurable_output', p);
            obj.state_x = [init_state.x; input.scheme_parameters.midrange(1,1)];
            obj.state_y = [init_state.y; input.scheme_parameters.midrange(2,1)];
            obj.w_bar = input.scheme_parameters.midrange;
            obj.max_increment(1,1) = 2 * input.scheme_parameters.alpha * input.scheme_parameters.dist_range(1,1) * (exp(eta * delta) - 1);
            obj.max_increment(2,1) = 2 * input.scheme_parameters.alpha * input.scheme_parameters.dist_range(2,1) * (exp(eta * delta) - 1);
            obj.is_saturator_active = true;


        end

        function obj = update(obj, x_measurements, y_measurements, control_input)
     
            obj.state_x = obj.A_upd_obs * obj.state_x + obj.B_upd_obs * control_input(1,1) + ...
                          obj.G' * (x_measurements - obj.C_measurable_output * obj.state_x);
            obj.state_y = obj.A_upd_obs * obj.state_y + obj.B_upd_obs * control_input(2,1) + ...
                          obj.G' * (y_measurements - obj.C_measurable_output * obj.state_y);

            if obj.is_saturator_active
                obj.w_bar(1,1) = quadprog(1, -obj.state_x(4), ...
                                        [1; -1; 1; -1], ...
                                        [obj.input.scheme_parameters.midrange(1,1) + obj.input.scheme_parameters.dist_range(1,1); ...
                                         - obj.input.scheme_parameters.midrange(1,1) + obj.input.scheme_parameters.dist_range(1,1); ...
                                         obj.w_bar(1,1) + obj.max_increment(1,1); ...
                                        -obj.w_bar(1,1) + obj.max_increment(1,1)]);
                obj.w_bar(2,1) = quadprog(1, -obj.state_y(4), ...
                                        [1; -1; 1; -1], ...
                                        [obj.input.scheme_parameters.midrange(2,1) + obj.input.scheme_parameters.dist_range(2,1); ...
                                         - obj.input.scheme_parameters.midrange(2,1) + obj.input.scheme_parameters.dist_range(2,1); ...
                                         obj.w_bar(2,1) + obj.max_increment(2,1); ...
                                        -obj.w_bar(2,1) + obj.max_increment(2,1)]);
            else 
                obj.w_bar(1,1) = obj.state_x(4,1);
                obj.w_bar(2,1) = obj.state_y(4,1);
            end

        end

        function w_bar = getDisturbance(obj)

                 w_bar = obj.w_bar;
                 
        end

    end

    properties (Access = private)

        A_upd_obs;
        B_upd_obs;
        G;
        C_measurable_output;
        state_x;
        state_y;
        w_bar;
        max_increment;
        input; 
        is_saturator_active;

    end

end