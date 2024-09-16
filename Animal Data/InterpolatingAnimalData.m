%% Neustadter 2007 Kinematic Data
% Last Updated: MJB 06/05/2024

clear all; close all; clc

folder = cd + "\Datasets\";

%% Getting Reference Length from High Res Image
dir = split(cd,"\");
dir = strjoin(dir(1:end-1),"\");
addpath(dir+"\Anatomical Parameter Initialization\")
mri_data = readcell("HighResolutionMidsagData.csv");

% scale bar
scale_bar_x = [mri_data{3:4,1}];
scale_bar_y = [mri_data{3:4,2}];
px_to_cm = sqrt(diff(scale_bar_x).^2 + diff(scale_bar_y).^2);

% radular stalk width
RSW_x = [mri_data{3:4,7}];
RSW_y = [mri_data{3:4,8}];
RSW_px = sqrt(diff(RSW_x).^2 + diff(RSW_y).^2);

RSW_per_mm = (px_to_cm / RSW_px) / 10;

%% Translation of the Odontophore
file = "Neustadter2007_Translation.csv";
data_trans = readmatrix(folder+file,'NumHeaderLines',2);

trans_bite_mean_data = data_trans(:,[1,2]);
trans_bite_p1std_data = data_trans(:,[3,4]);
trans_swallow_mean_data = data_trans(:,[5,6]);
trans_swallow_p1std_data = data_trans(:,[7,8]);

t_bite_peak_protract = data_trans(1,9);
t_swallow_peak_protract = data_trans(1,11);


[trans_bite_mean,trans_bite_std] = Kinematic_Fitting(trans_bite_mean_data,trans_bite_p1std_data);
[trans_swallow_mean,trans_swallow_std] = Kinematic_Fitting(trans_swallow_mean_data,trans_swallow_p1std_data);

%% Angle from Jaw Line to I6 Line
file = "Neustadter2007_JL_I6_Angle.csv";
data_JL_I6 = readmatrix(folder+file,'NumHeaderLines',2);

JL_I6_bite_mean_data = data_JL_I6(:,[1,2]);
JL_I6_bite_p1std_data = data_JL_I6(:,[3,4]);
JL_I6_swallow_mean_data = data_JL_I6(:,[5,6]);
JL_I6_swallow_p1std_data = data_JL_I6(:,[7,8]);

[JL_I6_bite_mean,JL_I6_bite_std] = Kinematic_Fitting(JL_I6_bite_mean_data,JL_I6_bite_p1std_data);
[JL_I6_swallow_mean,JL_I6_swallow_std] = Kinematic_Fitting(JL_I6_swallow_mean_data,JL_I6_swallow_p1std_data);

%% Angle from I6 to Radular Stalk Axis
file = "Neustadter2007_RS_I6_Angle.csv";
data_RS_I6 = readmatrix(folder+file,'NumHeaderLines',2);

RS_I6_bite_mean_data = data_RS_I6(:,[1,2]);
RS_I6_bite_p1std_data = data_RS_I6(:,[3,4]);
RS_I6_swallow_mean_data = data_RS_I6(:,[5,6]);
RS_I6_swallow_p1std_data = data_RS_I6(:,[7,8]);

[RS_I6_bite_mean,RS_I6_bite_std] = Kinematic_Fitting(RS_I6_bite_mean_data,RS_I6_bite_p1std_data);
[RS_I6_swallow_mean,RS_I6_swallow_std] = Kinematic_Fitting(RS_I6_swallow_mean_data,RS_I6_swallow_p1std_data);


%% Angle from Jaw Line to Radular Stalk Axis

RS_JL_bite_mean = @(t) (JL_I6_bite_mean(t) + RS_I6_bite_mean(t));
RS_JL_bite_std = @(t) sqrt(RS_I6_bite_std(t).^2 + JL_I6_bite_std(t).^2);
RS_JL_swallow_mean = @(t) (JL_I6_swallow_mean(t) + RS_I6_swallow_mean(t));
RS_JL_swallow_std = @(t) sqrt(RS_I6_swallow_std(t).^2 + JL_I6_swallow_std(t).^2);

%% Translation in Rejection (and calibration of RSW to mm)
file = "Novakovic2006_Translation.csv";
data_trans = readmatrix(folder+file,'NumHeaderLines',2);

trans_reject_mean_data = data_trans(:,[1,2]);
t_reject_peak_protract = data_trans(1,5) / max(trans_reject_mean_data(:,1));

trans_reject_mean_data(:,1) = trans_reject_mean_data(:,1) / max(trans_reject_mean_data(:,1));
trans_reject_p1std_data = trans_reject_mean_data; % no std data is presented
trans_swallow_mean_data_Nov = data_trans(:,[3,4]);
trans_swallow_mean_data_Nov(:,1) = trans_swallow_mean_data_Nov(:,1) / max(trans_swallow_mean_data_Nov(:,1));
trans_swallow_p1std_data_Nov = trans_swallow_mean_data_Nov;

[trans_reject_mean,~] = Kinematic_Fitting(trans_reject_mean_data,trans_reject_p1std_data);
[trans_swallow_mean_Nov,~] = Kinematic_Fitting(trans_swallow_mean_data_Nov,trans_swallow_p1std_data_Nov);

Nov_to_Neu = (max(trans_swallow_mean_data_Nov(:,2)) - min(trans_swallow_mean_data_Nov(:,2)))/(max(trans_swallow_mean_data(:,2)) - min(trans_swallow_mean_data(:,2)));

%% Rejection Timing Data
file = "Ye2006_RejectionTiming.csv";
data_reject_time = readmatrix("Datasets\"+file);
t_data = data_reject_time(:,1);

t_i2_on = t_data(1);                    % time I2 turns on (start of cycle)
t_i2_off = t_data(end-1);               % time I2 turns off
dur_i2 = t_i2_off - t_i2_on;            % measured duration of I2

t_rescaler = 3.0 / dur_i2;              % scalar to match graphical data to reported data

t_data = t_data * t_rescaler;           % adjust scaling of rejection data

t_i2_on = t_data(1);                    % time I2 turns on (start of cycle)
sig_i2_on = t_data(2) - t_i2_on;        % std on I2 turn on
t_i2_off = t_data(end-1);               % time I2 turns off
sig_i2_off = t_data(end) - t_i2_off;    % std on I2 turn off
t_end = t_data(3);                      % time cycle ends
sig_end = abs(t_data(4)-t_data(3));     % std on end of cycle

dur_i2 = t_i2_off - t_i2_on;
sig_dur_i2 = sqrt(sig_i2_on^2 + sig_i2_off^2);

reject_time_mean = t_end - t_i2_on;
reject_time_std = sqrt(sig_end^2 + sig_i2_on^2);

reject_perc_prot_mean = dur_i2 / reject_time_mean;
reject_perc_prot_std = reject_perc_prot_mean*sqrt((sig_dur_i2/dur_i2)^2 + (reject_time_std/reject_time_mean)^2);


%% Length of Swallowed Seaweed
Lum_SwallowingLength

%% Force on Seaweed
Gill_SwallowingForce

%% Animal Data Struct
t = linspace(0,1,1000);
animal_data = struct();

% offset to the peak protract -- used to align model with animal data
animal_data.t_bite_peak_protract = t_bite_peak_protract;
animal_data.t_swallow_peak_protract = t_swallow_peak_protract;
animal_data.t_reject_peak_protract = t_reject_peak_protract;

% translation data
animal_data.bite_trans_mean = @(t) trans_bite_mean(t) ;% * RSW_per_mm ;
animal_data.bite_trans_std = @(t) trans_bite_std(t) ;%* RSW_per_mm ;
animal_data.swallow_trans_mean = @(t) trans_swallow_mean(t) ;%* RSW_per_mm ;
animal_data.swallow_trans_std = @(t) trans_swallow_std(t) ;%* RSW_per_mm ;
animal_data.reject_trans_mean = @(t) (trans_reject_mean(t) / Nov_to_Neu) ;%* RSW_per_mm; 

% radular stalk to jawline data
animal_data.bite_angle_mean = RS_JL_bite_mean;
animal_data.bite_angle_std = RS_JL_bite_std;
animal_data.swallow_angle_mean = RS_JL_swallow_mean;
animal_data.swallow_angle_std = RS_JL_swallow_std;

% external behavioral data
animal_data.swallow_length_mean = @(t) length_swallow_mean(t) * 10 ; % converting from cm to mm
animal_data.swallow_length_std = @(t) length_swallow_std(t) * 10 ;
animal_data.swallow_force_mean = mean_force;
animal_data.swallow_force_std = std_force;

%% Timing Metrics from Bootstrap

bootstrapped_data = load("BootstrappedAnimalData.mat").animal_data;

animal_data.bite_time_mean = mean(bootstrapped_data.bite_time_data);
animal_data.bite_time_std = bootstrapped_data.bite_time_std;
animal_data.bite_prot_frac_mean = mean(bootstrapped_data.bite_prot_frac_data);
animal_data.bite_prot_frac_std =  bootstrapped_data.bite_prot_frac_std;

animal_data.unloaded_swallow_time_mean = mean(bootstrapped_data.uswallow_time_data);
animal_data.unloaded_swallow_time_std = bootstrapped_data.uswallow_time_std;
animal_data.unloaded_swallow_prot_frac_mean = mean(bootstrapped_data.uswallow_prot_frac_data);
animal_data.unloaded_swallow_prot_frac_std = bootstrapped_data.uswallow_prot_frac_std;


animal_data.lswallow_p_increase_time_mean = mean(bootstrapped_data.lswallow_cycle_p_increase_data);
animal_data.lswallow_p_increase_time_std = bootstrapped_data.lswallow_cycle_p_increase_std;
animal_data.lswallow_prot_frac_mean = mean(bootstrapped_data.lswallow_prot_frac_data);
animal_data.lswallow_prot_frac_std = bootstrapped_data.lswallow_prot_frac_std;

animal_data.reject_time_mean = reject_time_mean;
animal_data.reject_time_std = reject_time_std;
animal_data.reject_prot_frac_mean = reject_perc_prot_mean;
animal_data.reject_prot_frac_std = reject_perc_prot_std;



%% Range of Motion Estimates
% Biting Translation
[max_bite,ind_max_bite] = max(animal_data.bite_trans_mean(t));
[min_bite,ind_min_bite] = min(animal_data.bite_trans_mean(t));

animal_data.bite_trans_rom_mean = max_bite - min_bite;
animal_data.bite_trans_rom_std = sqrt(animal_data.bite_trans_std(t(ind_max_bite))^2 + animal_data.bite_trans_std(t(ind_min_bite))^2  );

% Swallowing Translation
[max_swallow,ind_max_swallow] = max(animal_data.swallow_trans_mean(t));
[min_swallow,ind_min_swallow] = min(animal_data.swallow_trans_mean(t));

animal_data.swallow_trans_rom_mean = max_swallow - min_swallow;
animal_data.swallow_trans_rom_std = sqrt(animal_data.swallow_trans_std(t(ind_max_swallow))^2 + animal_data.swallow_trans_std(t(ind_min_swallow))^2  );

% Rejectin Translation
[max_reject,ind_max_reject] = max(animal_data.reject_trans_mean(t));
[min_reject,ind_min_reject] = min(animal_data.reject_trans_mean(t));

animal_data.reject_trans_rom_mean = max_reject - min_reject;


% Biting Rotation
[max_bite,ind_max_bite] = max(animal_data.bite_angle_mean(t));
[min_bite,ind_min_bite] = min(animal_data.bite_angle_mean(t));

animal_data.bite_angle_rom_mean = max_bite - min_bite;
animal_data.bite_angle_rom_std = sqrt(animal_data.bite_angle_std(t(ind_max_bite))^2 + animal_data.bite_angle_std(t(ind_min_bite))^2  );

% Swallowing Translation
[max_swallow,ind_max_swallow] = max(animal_data.swallow_angle_mean(t));
[min_swallow,ind_min_swallow] = min(animal_data.swallow_angle_mean(t));

animal_data.swallow_angle_rom_mean = max_swallow - min_swallow;
animal_data.swallow_angle_rom_std = sqrt(animal_data.swallow_angle_std(t(ind_max_swallow))^2 + animal_data.swallow_angle_std(t(ind_min_swallow))^2  );


%% Saving Data
save("AnimalData","animal_data")
