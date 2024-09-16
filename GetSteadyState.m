%% GetSteadyState
%{
    Returns the steady state values of dynamic states and neurons of the
    given behavior.

%}

% old sensory variable names define the starting states of these states
% at time t_switch values switch from staring value X to 1-X


switch behavior
    case "biting"
        % Setting the stimuli parameters
        fixation_type = 0;
        sens_mechanical_lips = 1;
        sens_mechanical_grasper = 0;
        sens_chemical_lips = 1;

        CBI3_0 = 1;
        B31B32_0 = 1;
        B64_0 = 0;
        B6B9B3_0 = 0;

        t_switch_fixation = 1e6;
        t_switch_mech_lips = 1e6;
        t_switch_chem_lips = 1e6;
        t_switch_mech_grasp = 1e6;

    case "unloaded swallowing"
        % Setting the stimuli parameters
        fixation_type = 0;
        sens_mechanical_lips = 1;
        sens_mechanical_grasper = 1;
        sens_chemical_lips = 1;

        CBI3_0 = 1;
        B31B32_0 = 1;
        B64_0 = 0;
        B6B9B3_0 = 0;

        t_switch_fixation = 1e6;
        t_switch_mech_lips = 1e6;
        t_switch_chem_lips = 1e6;
        t_switch_mech_grasp = 1e6;
        
    case "loaded swallowing"
        % Setting the stimuli parameters
        fixation_type = 1;
        sens_mechanical_lips = 1;
        sens_mechanical_grasper = 1;
        sens_chemical_lips = 1;

        CBI3_0 = 1;
        B31B32_0 = 1;
        B64_0 = 0;
        B6B9B3_0 = 0;

        t_switch_fixation = 1e6;
        t_switch_mech_lips = 1e6;
        t_switch_chem_lips = 1e6;
        t_switch_mech_grasp = 1e6;
        
    case "rejection"
        % Setting the stimuli parameters
        fixation_type = 0;
        sens_mechanical_lips = 1;
        sens_mechanical_grasper = 1;
        sens_chemical_lips = 0;

        CBI3_0 = 0;
        B31B32_0 = 0;
        B64_0 = 1;
        B6B9B3_0 = 1;

        t_switch_fixation = 1e6;
        t_switch_mech_lips = 1e6;
        t_switch_chem_lips = 1e6;
        t_switch_mech_grasp = 1e6;
        
   otherwise
        error("ERROR: Behavior not available.")
end

% run the simulation and save the output struct
out = sim(simulinkFile); 

% find the rising edge of the B31/B32 motor pool to indicate start of cycle
dB31B32 = diff(out.B31B32); % getting steady state cycles
[~,starts] = findpeaks(dB31B32,"MinPeakDistance",1000); 
if length(starts) < 3
    i_ss = 1:length(dB31B32)+1;
else
    i_ss = starts(end-1):starts(end); % indices of steady state cycles
end

% starting time of the cycle
t0 = out.tout(i_ss(1));

% post-process to calculate the length of seaweed ingestion
CalculateLengthIngested;

