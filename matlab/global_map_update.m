function [new_map] = global_map_update(ARM_DATA_pos,old_map,global_damp_factor) %,move_history_angle_begin,move_history_angle_end,move_history_numsteps)
%
%Written by: Charlie Fisher, 2014-June-12
%
%This function performs motor map updating globally (the entire motor map)
%and locally (each region).
%
%-------------------------------------------------------------------------
%ARM_DATA_pos: This is an array of the positioners (in order) that are
%included in the data structure
%
%old_map: this is the old motor map (either slow or fast)
%
%global_damp_factor: an array that is the same size (and in same order as
%ARM_DATA_pos) that has the global damping factor for each motor.  These
%factors are generated from the move error vs move request plots (part of a
%different script)
%
%move_history_angle_begin: an array of the beginning move angle for every move
%made by the positioners
%
%move_history_angle_end: an array of the end move angle for every move
%made by the positioners
%
%move_history_numsteps: an array of the number of steps commanded for each
%move by the positioners
%-------------------------------------------------------------------------
%

%------Global Map Update-------------

num_Pos = length(ARM_DATA_pos);

for ii = 1:num_Pos
    ARM_ID = ['ARM_DATA_' num2str(ARM_DATA_pos(ii))];
    new_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes = old_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.*global_damp_factor.ARM_ID.Joint1_fwd;
    new_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes = old_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.*global_damp_factor.ARM_ID.Joint1_rev;
    new_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes = old_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.*global_damp_factor.ARM_ID.Joint2_fwd;
    new_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes = old_map.ARM_DATA.ARM_ID.SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.*global_damp_factor.ARM_ID.Joint2_rev;
end
%------------------------------------

%------Local Map Update--------------
