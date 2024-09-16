%% Code Controller for Aplysia Neuromechanical Model
% Last Updated: MB 09/15/2024
% 
% This script can be run to generate all of the results presented in 
% "Incorporating buccal mass planar mechanics and anatomical features 
% improves neuromechanical modeling of Aplysia feeding behavior."  Multiple 
% external scripts are called in order to keep this controller legible. 
% All external scripts, as well as the Simulink file, need to be in the 
% same folder so that the controller has access to them.

clear all; close all; clc

simulinkFile = "CMU_BORG_AplysiaModularNeuromechanics";
animal_data = load("Animal Data\AnimalData.mat").animal_data;

% Behaviors to Show -- if 1, the behavior will be simulated
show_biting = 1;
show_swallowing_unloaded = 1;
show_swallowing_loaded = 1;
show_rejection = 1;

% Comparisons and Saving
to_save = 1;                % indicator to specify if the simulation results should be saved (and overwrite the current data)
run_original_model = 1;     % indicator to specify if the Webster-Wood et al 2020 model should be run for comparison
recalc_data = true;        % indicator to specify whether or not to recalculate the damping convergence data (if false, loads saved data)

% plotting colors
present_model_color = [0.8,0.2,0.6];
boolean_model_color = [0.,.5,0.2];
animal_data_color = (0/255)*[1,1,1];
linewidth = 2;

%% Model Initialization

% this script calculates the necessary reference geometry, places model
% parameters in the workspace, and sets the initial values of the model
ReferenceGeometry_and_Initialization;

% vector of geometric parameters, in units of grasper radii (for lengths)
buccal_mass_params  =[1,L_head,0,H_lumen,xE6_anchor,theta_H,theta_fg];

MCC_0 = 1; % initializing metacerebral cell activity (1: model "hungry", 0: model "satiated") 
tend = 20; % setting end time of all simulations, ensures steady state reached

if run_original_model
    % get access to the Webster-Wood 2020 model code
    addpath(genpath('WebsterWood2020_Model'));
end

calculate_statistics = 1; % Do you want to calculate all of the statistics reported in the paper?

output_folder = 'ModelData\'; % location of saved data

%% Running New Model Simulations
seaweed_strength = 1e6; % seaweed is infinitely strong so the sim-slug cant break it

% Biting
if show_biting
    behavior = "biting";
    GetSteadyState; 
    out_bite = out;
    i_ss_bite = i_ss;
    t0b = t0;
    
    if to_save
        save(output_folder+"BiteOutput.mat","out_bite");
    end

end

% Unloaded Swallowing
if show_swallowing_unloaded
    behavior = "unloaded swallowing";
    GetSteadyState;
    out_uswallow = out;
    i_ss_uswallow = i_ss;
    t0us = t0;

    if to_save
        save(output_folder+"UnloadedSwallowOutput.mat","out_uswallow");
    end

end

% Loaded Swallowing 
if show_swallowing_loaded
    behavior = "loaded swallowing";
    GetSteadyState;
    out_lswallow = out;
    i_ss_lswallow = i_ss;
    t0ls = t0;

    if to_save
        save(output_folder+"LoadedSwallowOutput.mat","out_lswallow");
    end

end

% Rejection
if show_rejection
    behavior = "rejection";
    GetSteadyState;
    out_reject = out;
    i_ss_reject = i_ss;
    t0r = t0;

    if to_save
        save(output_folder+"RejectOutput.mat","out_reject");
    end

end

% Figure 5: Steady state simulation results
All_Behavior_Comparison

% Calculating kinematics data for Figure 7, Table 2, Figure 8, Table 3,
% Figure 9, and the statistical comparison
Calculating_Gross_Kinematics

%% Behavioral Switching

% Setting the stimuli parameters
fixation_type = 0;
sens_mechanical_lips = 1;
sens_mechanical_grasper = 0;
sens_chemical_lips = 1;

t_switch_fixation = 30;    % switch from unloaded to loaded swallowing at t=30 s
t_switch_mech_lips = 1e6;  % mechanical stimulation at the lips should never switch
t_switch_chem_lips = 54;   % switch from loaded swallowing to rejection at t=54 s
t_switch_mech_grasp = 10;  % switch to unloaded swallowing at t=10 s

tend = 80;                 % length of the simulation

% run the simulation
out = sim(simulinkFile); 
% post-processing to calculate length of ingested seaweed
CalculateLengthIngested;

out_behavior_switch = out;
if to_save
    save(output_folder+"BehaviorSwitching_Output.mat","out_behavior_switch");
end

% Figure 6: Behavioral switching
PlotFullBehaviorComparison

tend = 20; % resetting the length of simulations

%% Feeding on Breakable Seaweed

% Figure 7(a)
VaryingSeaweedStrength;

%% Running Webster-Wood Model Simulations

if run_original_model
    aplysia = AplysiaFeeding();

    if show_biting
        aplysia = aplysia.setSensoryStates('bite');
        aplysia_bite = aplysia.runSimulation();
        dB31B32 = diff(aplysia_bite.B31B32); % getting steady state cycles
        [~,starts] = findpeaks(dB31B32,"MinPeakDistance",100); 
        i_ss_bite_bool = starts(end-2):starts(end-1); % indices of steady state cycles
        tb_bool = aplysia_bite.StartingTime:aplysia_bite.TimeStep:aplysia_bite.EndTime;
        tb_bool = tb_bool';
        t0b_bool = tb_bool(i_ss_bite_bool(1));
    end

    if show_swallowing_unloaded
        aplysia = aplysia.setSensoryStates('swallow');
        aplysia_uswallow= aplysia.runSimulation();
        dB31B32 = diff(aplysia_uswallow.B31B32); % getting steady state cycles
        [~,starts] = findpeaks(dB31B32,"MinPeakDistance",100); 
        i_ss_uswallow_bool = starts(end-2):starts(end-1); % indices of steady state cycles
        tus_bool = aplysia_uswallow.StartingTime:aplysia_uswallow.TimeStep:aplysia_uswallow.EndTime;
        tus_bool = tus_bool';
        t0us_bool = tus_bool(i_ss_uswallow_bool(1));
    end

    if show_rejection
        aplysia = aplysia.setSensoryStates('reject');
        aplysia_reject = aplysia.runSimulation();
        dB31B32 = diff(aplysia_reject.B31B32); % getting steady state cycles
        [~,starts] = findpeaks(dB31B32,"MinPeakDistance",100); 
        i_ss_reject_bool = starts(end-2):starts(end-1); % indices of steady state cycles
        tr_bool = aplysia_reject.StartingTime:aplysia_reject.TimeStep:aplysia_reject.EndTime;
        tr_bool = tr_bool';
        t0r_bool = tr_bool(i_ss_reject_bool(1));
    end

end

%% Plotting Comparison with Neustadter 2007 Data

Neustadter2007_Comparison_NonNormalized
if run_original_model
    Neustadter2007_Comparison_Normalized
end

%% Plotting Comparison with Gill 2020 Data
Gill2020_Comparison

%% Plotting Comparison with Lum 2006 Data
Lum2006_Comparison

%% Damping Parameter Convergence Analysis
damping_data_file = "DampingParameterConvergenceData.mat";
DampingParameterConvergence

ReferenceGeometry_and_Initialization % run to reset parameters
MCC_0 = 1;
%% Equivalence and t Tests 

Statistical_Comparison

%% Plotting Frames of Model

PlotFrames

