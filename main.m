%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title: Robust Gait Generation Framework
% Author: Filippo M. Smaldone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clf; clear all;

%% input structure defining the gait task
input = struct;
input.delta = 0.01;
input.total_step_number = 10;
input.footstep_plan = struct;
input.footstep_plan.positions = zeros(input.total_step_number,3);
input.footstep_plan.orientations = zeros(input.total_step_number,3);
input.footstep_plan.timings = zeros(input.total_step_number,1);
input.footstep_plan.running_steps = zeros(input.total_step_number,1);
input.sim_time = 10;


%% state structure
state = struct;
state.x = zeros(3,1);
state.y = zeros(3,1);
state.support_foot = zeros(3,1);
state.base_orient = eye(3);


%% log structure
logs = struct;
logs.x_store = zeros(3, floor(input.sim_time/input.delta)); % x_c, x_dot_c, x_z
logs.y_store = zeros(3, floor(input.sim_time/input.delta)); % y_c, y_dot_c, y_z
logs.w_bar = zeros(2, floor(input.sim_time/input.delta)); % w_bar_x, w_bar_y
logs.w = zeros(2, floor(input.sim_time/input.delta)); % w_x, w_y
logs.actual_footsteps = zeros(3, input.total_step_number); % x, y, z


%% initialize the Robust Gait Generation Framework
first = GaitGeneration;

%% request back-up maneuver to initialize the controller (simulated)


%% simulation cycle
for iter = 1 : floor(input.sim_time/input.delta)
    
    % get measurement (simulate pertuebations)
    
    % solve step of gait generation algorithm
    
    % store logs
    
end



