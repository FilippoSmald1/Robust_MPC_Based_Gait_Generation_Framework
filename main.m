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


