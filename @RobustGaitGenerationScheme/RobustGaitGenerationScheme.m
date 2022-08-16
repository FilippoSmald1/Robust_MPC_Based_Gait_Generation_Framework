classdef RobustGaitGenerationScheme

    methods (Access = public)

        function obj = RobustGaitGenerationScheme(input, state)
            % constructor
            obj.a = 2;
            obj.test = 2;
            obj.input = input;
            obj.sm_instance = StandardMode(obj.input);
            obj.rm_instance = RecoveryMode(obj.input);
            obj.dob_instance = DisturbanceObserver(obj.input, state);  
        end

        function obj = update(obj, state)
             % control cycle
             zmpdot = zeros(2,1);
             obj.dob_instance.update([state.x(1,1); state.x(3,1)], ...
                                     [state.y(1,1); state.y(3,1)], ...
                                     zmpdot);            

        end

        function w_bar = getDisturbance(obj)

            w_bar = obj.dob_instance.getDisturbance();

        end

    end
    
    properties (Access = public)
        
        a = 1;
        test = 1;
        
    end

    properties (Access = private)
        
        sm_instance;
        rm_instance;
        dob_instance;
        input;
        
    end
    
end