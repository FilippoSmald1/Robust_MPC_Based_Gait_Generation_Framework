%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title: Robust Gait Generation Framework
% Author: Filippo M. Smaldone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clf; clear all; close all;


%% input data for the scheme
input = struct;

% parameters of the scheme
input.scheme_parameters = struct;
input.scheme_parameters.delta = 0.01; % sampling time
input.scheme_parameters.h = 0.75; % CoM height for LIP model
input.scheme_parameters.g = 9.81; % gravity acceleration constant
input.scheme_parameters.eta = sqrt(input.scheme_parameters.g / input.scheme_parameters.h);
input.scheme_parameters.T_c = 0.7; % prediction horizon
input.scheme_parameters.T_p = 2.5; % preview horizon
input.scheme_parameters.C = floor(input.scheme_parameters.T_c / input.scheme_parameters.delta);
input.scheme_parameters.P = floor(input.scheme_parameters.T_p / input.scheme_parameters.delta);
input.scheme_parameters.M = 2; % optimized footstep
input.scheme_parameters.F = 4; % available footstep from the plan at each time
input.scheme_parameters.midrange = [0.15; 0.15]; % (w_mx, w_my) in m/s^2
input.scheme_parameters.dist_range = [0.15; 0.15];  % (Deltaw_mx, Deltaw_my) in m/s^2
input.scheme_parameters.alpha = 0.25;
input.scheme_parameters.mi_max = 0.000025;
input.scheme_parameters.gamma_max = 0.00099;
input.scheme_parameters.epsilon = 0.005;
input.scheme_parameters.d_z = 0.1; % support polygon square width
input.scheme_parameters.d_ax = 0.5; % kinematic admissible region x dimension
input.scheme_parameters.d_ay = 0.25; % kinematic admissible region y dimension
input.scheme_parameters.ell_x = 0.0; % kinematic admissible region y displacement
input.scheme_parameters.ell_y = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.ell = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.d_ax_subsequent = 0.5; % kinematic admissible region x dimension
input.scheme_parameters.d_ay_subsequent = 0.25; % kinematic admissible region y dimension
input.scheme_parameters.ell_subsequent = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.ell_x_subsequent = 0.0; % kinematic admissible region y displacement
input.scheme_parameters.ell_y_subsequent = 0.2; % kinematic admissible region y displacement
input.scheme_parameters.footstep_weight_in_cost_function = 10; % 10 %1000000;
input.scheme_parameters.zmp_track_in_cost_function = 0.01;
input.scheme_parameters.v_max = 3;


%% handling non-convex constraints
% we begin with the standard convex KAR, plus two subregions that enable
% leg crossing maneuvers
input.kar = struct;
input.kar.number_of_subregions = 3;
input.kar.subregion_parameters = [input.scheme_parameters.d_ax, input.scheme_parameters.ell_x, input.scheme_parameters.d_ay, input.scheme_parameters.ell_y; ...
                                  0.15, 0.175, 0.15, 0.00; ...
                                  0.15, -0.175, 0.15, 0.00];


%% footstep plan
input.footstep_plan = struct;
input.footstep_plan.total_step_number = 18;
input.footstep_plan.positions = zeros(input.footstep_plan.total_step_number + 4,3);
input.footstep_plan.orientations = zeros(input.footstep_plan.total_step_number + 4,3);
input.footstep_plan.timings = zeros(input.footstep_plan.total_step_number + 4,1);
input.footstep_plan.running_steps = zeros(input.footstep_plan.total_step_number + 4,1);
input.footstep_plan.ds_duration = 0.2; % it is convenient to set a fixed duration for the double support
                                       % this can still be modified by the Step Timing Adaptation module
input.footstep_plan.ds_samples = floor(input.footstep_plan.ds_duration / input.scheme_parameters.delta);                                       
input.footstep_plan.starting_sf = "right";
input.footstep_plan.tail_x = zeros(input.scheme_parameters.P - input.scheme_parameters.C, 1);
input.footstep_plan.tail_y = zeros(input.scheme_parameters.P - input.scheme_parameters.C, 1);
input.footstep_plan.zmp_centerline_x = zeros(input.scheme_parameters.C, 1);
input.footstep_plan.zmp_centerline_y = zeros(input.scheme_parameters.C, 1);
input.footstep_plan.mapping_buffer = zeros(2 * input.scheme_parameters.P, input.scheme_parameters.M + 1);
input.sim_time = 10;

% build a simple footstep plan in the world frame
stride_length_x = 0.2;
lateral_displacement_y = 0.09;

number_of_virtual_steps = 4;

for i = 1 : input.footstep_plan.total_step_number + number_of_virtual_steps
     
   % footstep positions
   if i > 1
       input.footstep_plan.positions(i, 1) = (i-2) * stride_length_x;
   end
   if input.footstep_plan.starting_sf == "right"
       input.footstep_plan.positions(i, 2) = (- 1) ^ (i) * lateral_displacement_y;
   else
       input.footstep_plan.positions(i, 2) = (- 1) ^ (i - 1) * lateral_displacement_y;
   end
   input.footstep_plan.positions(i, 3) = 0;

   % footstep orientations (not considering orientation for simplicity)
     % keep all to zero

   % timings
   input.footstep_plan.timings(i, 1) = 1 *  (i - 1);

   % mapping to define the running steps (not considering flight phases here)
     % keep all to zero

end

% trick: model the initial double support as a 
% square centered between the feet
%input.footstep_plan.positions(1, 2) = 0;

% print the footstep plan
disp('input.footstep_plan.positions - (x,y,z) [m]')
disp(input.footstep_plan.positions)
disp('input.footstep_plan.timings [s]')
disp(input.footstep_plan.timings)



%% simulation parameters
simulation_parameters = struct;
simulation_parameters.delta = input.scheme_parameters.delta; %TODO: enable different operating frequencies
simulation_parameters.sim_time = 8;
simulation_parameters.sim_iter = 1;
simulation_parameters.sim_type = 'leg_crossing'; % 'basic_test', 'leg_crossing', 'obstacle'
simulation_parameters.obstacle_number = 1;
simulation_parameters.obstacles = [0.6 -0.09 0.05 0.05]; % x,y, size_x, size_y


%% state data 
state = struct;
state.x = zeros(3,1);
state.y = zeros(3,1);
state.w_bar = zeros(2,1);
state.sf_pos = input.footstep_plan.positions(1, 1:3)'; % position of the current support foot
state.sf_pos_ss  = input.footstep_plan.positions(1, 1:3)';
state.next_sf_pos = zeros(3,1);
state.next_sf_pos_ss = zeros(3,1);
state.current_sf = input.footstep_plan.starting_sf;
state.feasibility_region = [0; 0; 0; 0; 0; 0; 0; 0];
state.base_orient = eye(3);
state.footstep_counter = 1; % to query data from the plan
state.footstep_counter_rm = state.footstep_counter + 1;
state.footstep_counter_sm = state.footstep_counter;
state.step_time_iter = 1;
state.world_time_iter = 1;
state.sim_iter = 1;


%% log data
logs = struct;
logs.x_store = zeros(3, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % x_c, x_dot_c, x_z
logs.y_store = zeros(3, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % y_c, y_dot_c, y_z
logs.w_bar = zeros(2, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % w_bar_x, w_bar_y
logs.w = zeros(2, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % w_x, w_y
logs.actual_footsteps = zeros(3, input.footstep_plan.total_step_number); % x, y, z
logs.feasibility_region = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta));
logs.feasibility_region_2 = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta));
logs.feasibility_region_3 = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta));
logs.feasibility_region_full_logs = zeros(8, floor(simulation_parameters.sim_time/simulation_parameters.delta), 6);


%% initialize plotter
plotter = Plotter(logs, input, simulation_parameters);


%% CONTROLLER & SIMULATION START OPERATING HERE:
%
%
%
%
%% initialize the Robust Gait Generation Framework
wpg = RobustGaitGenerationScheme(input, state, simulation_parameters); % walking pattern generator


%% request back-up maneuver to initialize the controller (simulated)
% On a physical or multi-body simulated robot, this maneuver would consist in
% moving the CoM towards a feasible initialization for the MPC scheme:
% to do so, a simple polynomial interpolation from the current to the
% target CoM position is sufficient (to be tracked, e.g., via inverse kinematics).
% Here, we simply set the state to a proper initial value.
init_state_proposal = wpg.proposeFeasibleInitialState(state);
state.x(1,1) = init_state_proposal(1,1);
state.y(1,1) = init_state_proposal(3,1);
state.x(3,1) = state.sf_pos(1,1);
state.y(3,1) = state.sf_pos(2,1);


%% simulation cycle
for sim_iter = 1 : floor(simulation_parameters.sim_time / simulation_parameters.delta)

    % update iteration
    simulation_parameters.sim_iter = sim_iter;

    % get measurement (simulate pertuebations)
    state.x(2, 1) = state.x(2, 1) + input.scheme_parameters.delta * (0.15 + 0.1 * sin(2*pi*sim_iter*0.01/3)) ;
    state.y(2, 1) = state.y(2, 1) + input.scheme_parameters.delta * (0.15 + 0.1 * sin(2*pi*sim_iter*0.01/5)) ;    
    
    % disturbance and simple STA
    if strcmp(simulation_parameters.sim_type, 'basic_test')
        if sim_iter >= 230 && sim_iter < 240
            state.x(2, 1) = state.x(2, 1) + input.scheme_parameters.delta * 2.8 ;
            state.y(2, 1) = state.y(2, 1) - input.scheme_parameters.delta * 3.5 ;
        end
    end
    if strcmp(simulation_parameters.sim_type, 'leg_crossing')
        if sim_iter >= 230 && sim_iter < 240
            state.x(2, 1) = state.x(2, 1) + input.scheme_parameters.delta * 2.5 ;
            state.y(2, 1) = state.y(2, 1) + input.scheme_parameters.delta * 1.5 ;
        end
    end
    if strcmp(simulation_parameters.sim_type, 'obstacle')
        if sim_iter >= 330 && sim_iter < 340
            state.x(2, 1) = state.x(2, 1) + input.scheme_parameters.delta * 2 ;
            state.y(2, 1) = state.y(2, 1) + input.scheme_parameters.delta * 2 ;
        end
    end    
    % solve step of gait generation algorithm
    state = wpg.update(state);
    state.sim_iter = sim_iter;

    % store logs
    logs.x_store(:, sim_iter) = state.x;
    logs.y_store(:, sim_iter) = state.y;
    logs.w_bar(:, simulation_parameters.sim_iter) = wpg.getDisturbance();
    logs.actual_footsteps(:, state.footstep_counter) = state.sf_pos;
    logs.feasibility_region(:, sim_iter) = state.feasibility_region(:, 1);
    
    for k = 1 : size(state.feasibility_region, 2)
        logs.feasibility_region_full_logs(:, sim_iter, k) = state.feasibility_region(:, k);  
    end
    

    if sim_iter == 3 && false
        plotter.plotLogs(logs, state);
    end
    if mod(sim_iter, 50) == 0 && false
        plotter.plotLogsAtTimeK(logs, state, sim_iter);
    end
    
    if state.footstep_counter >= 3
            1;    
    end  
    
    str = num2str( sim_iter / floor(simulation_parameters.sim_time / simulation_parameters.delta) * 100 );
    str = strcat('simulation is at', " ", str, '%');
    disp(str)
    
end
%
%
%
%
%% CONTROLLER & SIMULATION STOP OPERATING HERE

%% plot the logs at a certain time 
t_k = 6.5;
plotter.plotLogsAtTimeK(logs, state, floor(t_k / simulation_parameters.delta));

