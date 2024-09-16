%% ~~~~~~~~~~~~~~~ Calculating Gross Kinematics ~~~~~~~~~~~~~~~ 
%{
    This script calculates all of the gross kinematic values that are used
    for statistical comparison with the animal data.

    Variable Naming Convention:
    - t_cycle_{behavior}: [s] length of steady state cycle
    - {behavior}_percent_protraction: [%] % of cycle spent in protraction (ie the odontophore is moving forward)
    - {behavior}_percent_retraction: [%] % of cycle spent in retraction (ie the odontophore is moving backwards)
    - {behavior}_x_ROM_mm: [mm] translational range of motion
    - {behavior}_theta_ROM: [deg] rotational range of motion
%}

model_data = struct();

%% Biting

% Rotation
t_bite_ss = out_bite.tout(i_ss_bite)-t0b; % steady state time vector
I2_bite = t_bite_ss(out_bite.B31B32(i_ss_bite) > 0);
t_cycle_bite = t_bite_ss(end); % [s] length of steady state cycle

bite_percent_protraction = I2_bite(end) / t_cycle_bite; 
bite_percent_retraction = 1 - bite_percent_protraction; 

% translational range of motion
bite_x_RoM = max(out_bite.x_gh(i_ss_bite)) - min(out_bite.x_gh(i_ss_bite));
bite_x_RoM_mm = bite_x_RoM*R_ref*10*(0.5*(L0_I1v + L0_I1d)); % R_ref brings back real units, 10 scales cm to mm, 0.5*(L0_I1v+L0_I1d) undoes normalization

% angular range of motion
bite_theta_RoM = (180/pi)*( max(out_bite.theta_g_animal(i_ss_bite)) - min(out_bite.theta_g_animal(i_ss_bite)) );

% Storing Values in struct
model_data.bite_time = t_cycle_bite; % [s]
model_data.bite_prot_frac = bite_percent_protraction; % [s/s]
model_data.bite_ret_frac = bite_percent_retraction; % [s/s]
model_data.bite_trans_rom = bite_x_RoM_mm;  % [mm]
model_data.bite_angle_rom = bite_theta_RoM; % [deg]
%% Unloaded Swallowing

% Rotation
t_uswallow_ss = out_uswallow.tout(i_ss_uswallow)-t0us;
I2_uswallow =  t_uswallow_ss(out_uswallow.B31B32(i_ss_uswallow) > 0);
t_cycle_uswallow = t_uswallow_ss(end);

uswallow_percent_protraction = I2_uswallow(end) / t_cycle_uswallow; 
uswallow_percent_retraction = 1 - uswallow_percent_protraction; 

uswallow_x_RoM = max(out_uswallow.x_gh(i_ss_uswallow)) - min(out_uswallow.x_gh(i_ss_uswallow));
uswallow_x_RoM_mm = uswallow_x_RoM*R_ref*10*(0.5*(L0_I1v + L0_I1d)); % R_ref brings back real units, 10 scales cm to mm, 0.5*(L0_I1v+L0_I1d) undoes normalization

uswallow_theta_RoM = (180/pi)*( max(out_uswallow.theta_g_animal(i_ss_uswallow)) - min(out_uswallow.theta_g_animal(i_ss_uswallow)) );

% Storing Values in struct
model_data.uswallow_time = t_cycle_uswallow; % [s]
model_data.uswallow_prot_frac = uswallow_percent_protraction; % [s/s]
model_data.uswallow_ret_frac = uswallow_percent_retraction; % [s/s]
model_data.uswallow_trans_rom = uswallow_x_RoM_mm;  % [mm]
model_data.uswallow_angle_rom = uswallow_theta_RoM; % [deg]
%% Loaded Swallowing

% Rotation
t_lswallow_ss = out_lswallow.tout(i_ss_lswallow)-t0ls;
I2_lswallow =  t_lswallow_ss(out_lswallow.B31B32(i_ss_lswallow) > 0);
t_cycle_lswallow = t_lswallow_ss(end);

lswallow_percent_protraction = I2_lswallow(end) / t_cycle_lswallow;
lswallow_percent_retraction = 1 - lswallow_percent_protraction;

lswallow_x_RoM = max(out_lswallow.x_gh(i_ss_lswallow)) - min(out_lswallow.x_gh(i_ss_lswallow));
lswallow_x_RoM_mm = lswallow_x_RoM*R_ref*10*(0.5*(L0_I1v + L0_I1d)); % R_ref brings back real units, 10 scales cm to mm, 0.5*(L0_I1v+L0_I1d) undoes normalization

lswallow_theta_RoM = (180/pi)*(max(out_lswallow.theta_g_animal(i_ss_lswallow)) - min(out_lswallow.theta_g_animal(i_ss_lswallow)));

% Storing Values in struct
model_data.lswallow_time = t_cycle_lswallow; % [s]
model_data.lswallow_prot_frac = lswallow_percent_protraction; % [s/s]
model_data.lswallow_ret_frac = lswallow_percent_retraction; % [s/s]
model_data.lswallow_trans_rom = lswallow_x_RoM_mm;  % [mm]
model_data.lswallow_angle_rom = lswallow_theta_RoM; % [deg]
%% Rejection

% Rotation
t_reject_ss = out_reject.tout(i_ss_reject)-t0r;
I2_reject =  t_reject_ss(out_reject.B31B32(i_ss_reject) > 0);
t_cycle_reject = t_reject_ss(end);

reject_percent_protraction = I2_reject(end) / t_cycle_reject; 
reject_percent_retraction = 1 - reject_percent_protraction;

reject_x_RoM = max(out_reject.x_gh(i_ss_reject)) - min(out_reject.x_gh(i_ss_reject));
reject_x_RoM_mm = reject_x_RoM*R_ref*10*(0.5*(L0_I1v + L0_I1d)); % R_ref brings back real units, 10 scales cm to mm, 0.5*(L0_I1v+L0_I1d) undoes normalization

reject_theta_RoM = (180/pi)*( max(out_reject.theta_g_animal(i_ss_reject)) - min(out_reject.theta_g_animal(i_ss_reject)) );

% Storing Values in struct
model_data.reject_time = t_cycle_reject; % [s]
model_data.reject_prot_frac = reject_percent_protraction; % [s/s]
model_data.reject_ret_frac = reject_percent_retraction; % [s/s]
model_data.reject_trans_rom = reject_x_RoM_mm;  % [mm]
model_data.reject_angle_rom = reject_theta_RoM; % [deg]


%% Saving struct for statistical analysis

if to_save
    save("Model_SummaryMetrics.mat","model_data")
end
