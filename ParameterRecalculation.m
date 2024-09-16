%% Parameter Recalculation with New Time Scalar
%{
    Code is copied directly from ReferenceGeometry_and_Initialization.m
    starting from line 127
%}

% damping parameters
C_yg = C_xg;
C_theta = C_xg/(2*pi);
C_xd = C_xg;
C_xv = C_xd;
C_xh = C_xg;

% frictional parameters
usg = 1; % from original model
ukg = 0.95*usg;
ush = 1;
ukh = 0.75*ush;

% anatomical parameters
x_gh_ref = xg0 - xh0;
x_gh_max = (L_head - grasperR - x_gh_ref);

% MUSCLE PARAMETERS

% I1d
AI1d0 = 0.0;            % [ ] starting activation
TI1d0 = 0;              % [ ] starting normalized tension
tau_I1d = time_scaler*time_constants.tau_I3;          % [1/s] muscle time constant
L0_I1d = L1d0;          % [R] rest length of I1d
Lmax_I1d = L1d0*(1.45);  % [R] max length of I1d
I1d_params = [AI1d0,TI1d0,tau_I1d,L0_I1d,Lmax_I1d]; % muscle params of I1d
FI1d_max = force_scaler*3.2;

% I1v
AI1v0 = 0.0;            % [ ] starting activation
TI1v0 = 0;              % [ ] starting normalized tension
tau_I1v = tau_I1d;          % [1/s] muscle time constant
L0_I1v = L1v0;          % [R] rest length of I1v
Lmax_I1v = L1v0*(1.45);  % [R] max length of I1v
I1v_params = [AI1v0,TI1v0,tau_I1v,L0_I1v,Lmax_I1v]; % muscle params of I1v
FI1v_max = FI1d_max;

% I2 (for linear spring + 1st order, need to update for Hill)
AI20 = 0.;             % [ ] starting activation
TI20 = 0;               % [ ] starting normalized tension
tau_I2_ingestive = time_scaler*time_constants.tau_I2_ingestive; % [1/s] ingestive muscle time constant
tau_I2_egestive = time_constants.tau_I2_egestive;  % [1/s] egestive muscle time constant
L0_I2 = LI2*1;            % [R] rest length of I2
Lmax_I2 = L0_I2*(1.6);    % [R] max length of I2
I2_params = [AI20,TI20,tau_I2_ingestive,tau_I2_egestive,L0_I2,Lmax_I2]; % muscle params of I2
FI2_max = force_scaler*1;

% I3 
AI30 = 0.0;            % [ ] starting activation
TI30 = 0;              % [ ] starting normalized tension
tau_I3 = tau_I1v;          % [1/s] muscle time constant
I3_params = [AI30,TI30,tau_I3]; % muscle params of I3
FI3_max = force_scaler*3;

% I3 anterior
AI3_ant0 = 0.0;            % [ ] starting activation
TI3_ant0 = 0;              % [ ] starting normalized tension
tau_I3ant = time_scaler*time_constants.tau_I3ant;     % [1/s] muscle time constant
I3ant_params = [AI3_ant0,TI3_ant0,tau_I3ant]; % muscle params of I3 ant
FI3_ant_max = force_scaler*0.1;

% I4
AI40 = 0.0;            % [ ] starting activation
TI40 = 0;              % [ ] starting normalized tension
tau_I4 = time_scaler*time_constants.tau_I4;    % [1/s] muscle time constant
I4_params = [AI40,TI40,tau_I4]; % muscle params of I4
FI4_max = force_scaler*0.9;

% Hinge
A_hinge0 = 0.0;            % [ ] starting activation
T_hinge0 = 0;              % [ ] starting normalized tension
tau_hinge = time_scaler*time_constants.tau_hinge;          % [1/s] muscle time constant
L0_hinge = 0;%L_hinge;           % [R] rest length of hinge
Lmax_hinge = (7/5)*(L0_I1v - 0.5*grasperR); %8/15; %L_hinge*(1.4);   % [R] max length of hinge
Hinge_params = [A_hinge0,T_hinge0,tau_hinge,L0_hinge,Lmax_hinge]; % muscle params of hinge
FHinge_max = force_scaler*0.292;

% E1
AE10 = 0.0;            % [ ] starting activation
TE10 = 0;              % [ ] starting normalized tension
tau_E1 = time_scaler*0.606;          % [1/s] muscle time constant
L0_E1 = LE1;          % [R] rest length of E1
Lmax_E1 = LE1*(1.5);  % [R] max length of E1
E1_params = [AE10,TE10,tau_E1,L0_E1,Lmax_E1]; % muscle params of E1
FE1_max = force_scaler*0.2;

% E2
AE20 = 0.0;            % [ ] starting activation
TE20 = 0;              % [ ] starting normalized tension
tau_E2 = tau_E1;          % [1/s] muscle time constant
L0_E2 = LE2;          % [R] rest length of E2
Lmax_E2 = LE2*( Lmax_E1 / L0_E1 );  % [R] max length of E2
E2_params = [AE20,TE20,tau_E2,L0_E2,Lmax_E2]; % muscle params of E2
FE2_max = FE1_max;

% E6
AE60 = 0.0;            % [ ] starting activation
TE60 = 0;              % [ ] starting normalized tension
tau_E6 = tau_E1;          % [1/s] muscle time constant
L0_E6 = LE6;          % [R] rest length of E6
Lmax_E6 = LE6*( Lmax_E1 / L0_E1 );  % [R] max length of E6
E6_params = [AE60,TE60,tau_E6,L0_E6,Lmax_E6]; % muscle params of E6
FE6_max = FE1_max;

% Passive Springs
K_Eso = force_scaler*0.01;      % [N/R] esophagus spring stiffness
K_H = force_scaler*2;           % [N/R] head spring stiffness   
K_const = force_scaler*0.5;     % [N/R] penalty stiffness

%% Neuron Initial Conditions
CBI3_stimOFF_0 = 0;
CBI3_refractory_0 = 0;
B40B30_0 = 0;
MCC_0 = 0;
CBI2_0 = 0;
CBI3_stimON_0 = 0;
CBI4_0 = 0;
B64_0 = 0;
B4B5_0 = 0;
B20_0 = 0;
B40B30_offTime_0 = 0;
B8_0 = 0;
B38_0 = 0;
B6B9B3_0 = 0;
CBI3_0 = 1;
B31B32_0 = 1;
B7_0 = 0;
B43B45_0 = 0;
C1_0 = 0;
C2_0 = 0;
C6_0 = 0;


%% Neuron Model Parameters
params_struct = {};

% Biting Parameters
params_struct.thresh_B31_bite_off =            0.54;    % retraction level at which to turn I2 back on when its of
params_struct.thresh_B64_bite_protract =       0.56;    % protraction level at which to turn I3 on (once I2 turns off)
params_struct.thresh_B31_bite_on =             1.1 ;    % protraction level at which to turn I2 off
params_struct.thresh_B7_bite_protract =        0.80 ;   % protraction level at which to turn the hinge on
params_struct.thresh_B7_bite_pressure =        0.7;     % pressure level at which B7 turns on
params_struct.thresh_B6B9B3_bite_pressure =    0.5;     % pressure threshold needed before I3 turns on 
params_struct.thresh_B31_pressure_bite =       0.15;    % pressure threshold at which I2 can turn on

% Swallowing Parameters
params_struct.thresh_B64_swallow_protract =    0.24;    % protraction level above which to turn on I3
params_struct.thresh_B31_swallow_off =         0.26;    % retraction level at which to turn I2 back on 
params_struct.thresh_B31_swallow_on =          0.74;    % protraction level at which to turn I2 off
params_struct.thresh_B6B9B3_swallow_pressure = 0.75;    % pressure threshold needed before I3 turns on 
params_struct.thresh_B38_retract =             0.5;     % threshold below which B38 turns on
params_struct.thresh_B31_pressure_swallow =    0.25;    % pressure threshold at which I2 can turn on

% Rejection Parameters
params_struct.thresh_B64_reject_protract =     0.22;    % protraction level above which to turn on I3
params_struct.thresh_B4B5_protract =           1.10;    % protraction level above which to inhibit the I3 via high stimulation of B4/B5
params_struct.thresh_B31_reject_off =          0.4;     % retraction level at which to turn I2 back on when its of
params_struct.thresh_B31_reject_on =           0.81;    % protraction level at which to turn I2 off
params_struct.thresh_B7_reject_protract =      0.95;    % protraction level at which to turn the hinge on
params_struct.thresh_B7_reject_pressure =      0.6;     % pressure level at which B7 turns on
params_struct.thresh_B6B9B3_reject_pressure =  0.3;     % pressure threshold needed before I3 turns on 
params_struct.thresh_B31_pressure_reject =     0.09;    % pressure threshold at which I2 can turn on

params_struct.refractory_CBI3 = 6.06 * time_scaler;             % [s] duration of the CBI-3 refractory period
params_struct.postActivityExcitation_B40B30 = 3.8*time_scaler;  % [s] duration of B8 stimulation after inactivation of B40B30

% converting params to a vector for use in parallel
fields = fieldnames(params_struct);
param_vec = zeros(length(fields),1);
for i=1:length(fields)
    param_vec(i) = params_struct.(fields{i});
end

params = param_vec; 

