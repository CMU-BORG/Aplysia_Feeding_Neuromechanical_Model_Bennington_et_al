%% Fitting Time Constant Parameters from animal data

% clear all; close all; clc;
dir = split(cd,"\");
dir = strjoin(dir(1:end-1),"\");
addpath(dir+"\Animal Data\Datasets\");

ComparisonWithYu; % getting I2 time constant for double first order filter

params_from_data = struct();
tau_I2 = median(tau_fit); % I2 time constant to normalize

params_from_data.tau_I2_ingestive = 1; % normalized time constant
params_from_data.tau_I2_egestive = params_from_data.tau_I2_ingestive / median(beta_fit);

%% Bulk I3 and I1 time constant from average time constants in Sukhnandan 2024
tau_I3_mean = mean([0.6,0.6/0.15]); % time constant for nonlinear first order filter
tau_I3_DFO = tau_I3_mean * (tau_I2 / 2.45); % correction factor for nonlinear first order to DFO filter
params_from_data.tau_I3 = tau_I3_DFO / tau_I2; 


%% I3 pinch data - from McManus et al. 2014
McManusI3aDataFitting;
params_from_data.tau_I3ant = tau / tau_I2;
fprintf("    tau [normalized] = %.3f\n",tau / tau_I2)
saveas(gcf,"Figures\McManus2014.svg",'svg')

%% I4 data - from Morton et al. 1993
MortonI4DataFitting;
params_from_data.tau_I4 = tau / tau_I2;
fprintf("    tau [normalized] = %.3f\n",tau / tau_I2)
saveas(gcf,"Figures\Morton1993.svg",'svg')

%% Hinge Data - from Sutton et al 2004
SuttonHingeDataFitting;
params_from_data.tau_hinge = tau / tau_I2;
fprintf("    tau [normalized] = %.3f\n",tau / tau_I2)
saveas(gcf,"Figures\Sutton2004.svg",'svg')

%% Saving Data
save("AnimalTimeConstants","params_from_data")
