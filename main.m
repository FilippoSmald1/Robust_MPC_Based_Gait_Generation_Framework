%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Title: Robust Gait Generation Framework
% Author: Filippo M. Smaldone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clf; clear all;
addpath modes_of_operation;
addpath ancillary_classes;

%% input structure defining the gait task
input = struct;

% parameters of the scheme
input.scheme_parameters = struct;
input.scheme_parameters.delta = 0.01; % sampling time
input.scheme_parameters.h = 0.75; % CoM height for LIP model
input.scheme_parameters.g = 9.81; % gravity acceleration constant

% footstep plan
input.footstep_plan = struct;
input.footstep_plan.total_step_number = 10;
input.footstep_plan.positions = zeros(input.footstep_plan.total_step_number,3);
input.footstep_plan.orientations = zeros(input.footstep_plan.total_step_number,3);
input.footstep_plan.timings = zeros(input.footstep_plan.total_step_number,1);
input.footstep_plan.running_steps = zeros(input.footstep_plan.total_step_number,1);
input.sim_time = 10;

% build a simple footstep plan in the world frame
stride_length_x = 0.2;
lateral_displacement_y = 0.09;

for i = 1:input.footstep_plan.total_step_number 
    
   % footstep positions
   input.footstep_plan.positions(i, 1) = (i-1) * stride_length_x;
   input.footstep_plan.positions(i, 2) = (- 1) ^ (i - 1) * lateral_displacement_y;
   input.footstep_plan.positions(i, 3) = 0;
   
   % footstep orientations (TODO)
     % zeros
     
   % timings
   input.footstep_plan.timings(i, 1) = 0.7 * (i - 1);
   
   % mapping to define the running steps (for further developments)
     % zeros
   
end

% print the footstep plan
disp('input.footstep_plan.positions - (x,y,z) [m]')
disp(input.footstep_plan.positions)
disp('input.footstep_plan.timings [s]')
disp(input.footstep_plan.timings)
pause


%% simulation parameters
simulation_parameters = struct;
simulation_parameters.delta = input.scheme_parameters.delta; %TODO: enable different operating frequencies
simulation_parameters.sim_time = 5;


%% state structure
state = struct;
state.x = zeros(3,1);
state.y = zeros(3,1);
state.support_foot = zeros(3,1);
state.base_orient = eye(3);


%% log structure
logs = struct;
logs.x_store = zeros(3, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % x_c, x_dot_c, x_z
logs.y_store = zeros(3, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % y_c, y_dot_c, y_z
logs.w_bar = zeros(2, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % w_bar_x, w_bar_y
logs.w = zeros(2, floor(simulation_parameters.sim_time/simulation_parameters.delta)); % w_x, w_y
logs.actual_footsteps = zeros(3, input.footstep_plan.total_step_number); % x, y, z


%% initialize the Robust Gait Generation Framework
wpg = RobustGaitGenerationScheme; % walking pattern generator



%% CONTROLLER & SIMULATION START OPERATING HERE:
%
%
%
%
%% request back-up maneuver to initialize the controller (simulated)


%% simulation cycle
for iter = 1 : floor(simulation_parameters.sim_time / simulation_parameters.delta)
    
    % get measurement (simulate pertuebations)
    
    % solve step of gait generation algorithm
    
    % store logs
    
end
%
%
%
%
%% CONTROLLER & SIMULATION STOP OPERATING HERE

%% plot the logs


