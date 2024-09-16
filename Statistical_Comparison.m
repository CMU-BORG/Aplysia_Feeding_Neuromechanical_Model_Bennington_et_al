%% Statistical Comparison for Bootstrapped Animal Data
%{
    Performs both t-test and equivalence tests for the model simulation
    data compared to the aggregated animal data. For complete details, see
    Supporting Information 3.
%}


% Setting for the statistical tests

% equivalence bound in terms of percent of the mean, chosen to achieve a
% power of 0.8 for observed difference of 0. This number is based on the
% number of animals available for each behavior
d_bite = 1.27;      % Cullins, N = 7
d_uswallow = 0.92;  % Cullins + Gill, N = 12
d_lswallow = 1.65;  % Gill, N = 5
d_reject = 1.65;    % Ye, N = 5

bootstrap_animal_data = load("Animal Data\BootstrappedAnimalData.mat").animal_data;
model_data = load("Model_SummaryMetrics.mat").model_data;

figure("Position",[100,100,1000,500],"Color","w"); hold all

xline(0,'--k')
metric_ind = 7;
xlim(2*[-1,1])
labels = fliplr({"Biting Cycle Duration", ...
                 "Biting Percent Protraction", ...
                 "Unloaded Swallowing Cycle Duration", ...
                 "Unloaded Swallowing Percent Protraction", ...
                 "Loaded Cycle Duration Percent Increase", ...
                 "Loaded Cycle Percent Protraction", ...
                 "Rejection Cycle Duration", ...
                 "Rejection Percent Protraction"});

col_95 = [1,0,0];

%% I. Biting
fprintf("Biting Metric Comparison:\n")

%% I. a) Cycle Duration

model_bite_time = model_data.bite_time;
animal_bite_time_means = bootstrap_animal_data.bite_time_data;
animal_bite_time_std = bootstrap_animal_data.bite_time_std;
bite_time_diff = (animal_bite_time_means - model_bite_time);
bite_time_90CI = quantile(bite_time_diff,[0.05,0.95]);
bite_time_95CI = quantile(bite_time_diff,[0.025,0.975]);

plot(d_bite*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot(mean(bite_time_diff) / animal_bite_time_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(bite_time_95CI / animal_bite_time_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(bite_time_90CI / animal_bite_time_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

fprintf("\tCycle Time:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",mean(animal_bite_time_means),std(animal_bite_time_means),animal_bite_time_std)
fprintf("\t\tModel Value: %.3f\n",model_bite_time)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(mean(animal_bite_time_means) - model_bite_time)/mean(animal_bite_time_means),...
                                                        (mean(animal_bite_time_means) - model_bite_time)/animal_bite_time_std )

metric_ind = metric_ind - 1; 
%% I. b) Percent Protraction

model_bite_perc_prot = model_data.bite_prot_frac;
animal_bite_perc_prot_means = bootstrap_animal_data.bite_prot_frac_data;
animal_bite_perc_prot_std = bootstrap_animal_data.bite_prot_frac_std;
bite_perc_proc_frac_diff = (animal_bite_perc_prot_means - model_bite_perc_prot);
bite_perc_proc_90CI = quantile(bite_perc_proc_frac_diff,[0.05,0.95]);
bite_perc_proc_95CI = quantile(bite_perc_proc_frac_diff,[0.025,0.975]);

plot(d_bite*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot(mean(bite_perc_proc_frac_diff)/animal_bite_perc_prot_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(bite_perc_proc_95CI/animal_bite_perc_prot_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(bite_perc_proc_90CI/animal_bite_perc_prot_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

fprintf("\tPercent Protraction:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",mean(animal_bite_perc_prot_means),std(animal_bite_perc_prot_means),animal_bite_perc_prot_std)
fprintf("\t\tModel Value: %.3f\n",model_bite_perc_prot)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(mean(animal_bite_perc_prot_means) - model_bite_perc_prot)/mean(animal_bite_perc_prot_means),...
                                                        (mean(animal_bite_perc_prot_means) - model_bite_perc_prot)/ animal_bite_perc_prot_std)


metric_ind = metric_ind - 1;
%% II. Unloaded Swallowing
fprintf("Unloaded Swallowing Metric Comparison:\n")

%% I. a) Cycle Duration

model_uswallow_time = model_data.uswallow_time;
animal_uswallow_time_means = bootstrap_animal_data.uswallow_time_data;
animal_uswallow_time_std = bootstrap_animal_data.uswallow_time_std;
uswallow_time_diff = animal_uswallow_time_means - model_uswallow_time;
uswallow_perc_proc_90CI = quantile(uswallow_time_diff,[0.05,0.95]);
uswallow_perc_proc_95CI = quantile(uswallow_time_diff,[0.025,0.975]);

plot(d_uswallow*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot(mean(uswallow_time_diff)/animal_uswallow_time_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(uswallow_perc_proc_95CI/animal_uswallow_time_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(uswallow_perc_proc_90CI/animal_uswallow_time_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

fprintf("\tCycle Time:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",mean(animal_uswallow_time_means),std(animal_uswallow_time_means),animal_uswallow_time_std)
fprintf("\t\tModel Value: %.3f\n",model_uswallow_time)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(mean(animal_uswallow_time_means) - model_uswallow_time)/mean(animal_uswallow_time_means),...
                                                    (mean(animal_uswallow_time_means) - model_uswallow_time)/animal_uswallow_time_std)


metric_ind = metric_ind - 1;

%% II. b) Percent Protraction

model_uswallow_perc_prot = model_data.uswallow_prot_frac;
animal_uswallow_perc_prot_means = bootstrap_animal_data.uswallow_prot_frac_data;
animal_uswallow_time_std = bootstrap_animal_data.uswallow_prot_frac_std;
uswallow_perc_proc_frac_diff = (animal_uswallow_perc_prot_means - model_uswallow_perc_prot);
uswallow_perc_proc_90CI = quantile(uswallow_perc_proc_frac_diff,[0.05,0.95]);
uswallow_perc_proc_95CI = quantile(uswallow_perc_proc_frac_diff,[0.025,0.975]);

plot(d_uswallow*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot(mean(uswallow_perc_proc_frac_diff)/animal_uswallow_time_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(uswallow_perc_proc_95CI/animal_uswallow_time_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(uswallow_perc_proc_90CI/animal_uswallow_time_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

fprintf("\tPercent Protraction:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",mean(animal_uswallow_perc_prot_means),std(animal_uswallow_perc_prot_means),animal_uswallow_time_std)
fprintf("\t\tModel Value: %.3f\n",model_uswallow_perc_prot)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(mean(animal_uswallow_perc_prot_means) - model_uswallow_perc_prot)/mean(animal_uswallow_perc_prot_means),...
                                                        (mean(animal_uswallow_perc_prot_means) - model_uswallow_perc_prot)/animal_uswallow_time_std)


metric_ind = metric_ind - 1;

%% III. Loaded Swallows
fprintf("Loaded Swallowing Metric Comparison:\n")

%% III. a) Percent Increase in Cycle Duration

model_lswallow_cycle_p_increase = (model_data.lswallow_time - model_uswallow_time)/model_uswallow_time;
animal_lswallow_cycle_p_increase_means = bootstrap_animal_data.lswallow_cycle_p_increase_data;
animal_lswallow_cycle_p_increase_std = bootstrap_animal_data.lswallow_cycle_p_increase_std;
lswallow_cycle_p_increase_diff = animal_lswallow_cycle_p_increase_means - model_lswallow_cycle_p_increase;
lswallow_cycle_p_increase_90CI = quantile(lswallow_cycle_p_increase_diff,[0.05,0.95]);
lswallow_cycle_p_increase_95CI = quantile(lswallow_cycle_p_increase_diff,[0.025,0.975]);

plot(d_lswallow*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot(mean(lswallow_cycle_p_increase_diff)/animal_lswallow_cycle_p_increase_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(lswallow_cycle_p_increase_95CI/animal_lswallow_cycle_p_increase_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(lswallow_cycle_p_increase_90CI/animal_lswallow_cycle_p_increase_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

fprintf("\tPercent Increase in Cycle Duration:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",mean(animal_lswallow_cycle_p_increase_means),std(animal_lswallow_cycle_p_increase_means),animal_lswallow_cycle_p_increase_std)
fprintf("\t\tModel Value: %.3f\n",model_lswallow_cycle_p_increase)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(mean(animal_lswallow_cycle_p_increase_means) - model_lswallow_cycle_p_increase)/mean(animal_lswallow_cycle_p_increase_means),...
                                                        (mean(animal_lswallow_cycle_p_increase_means) - model_lswallow_cycle_p_increase)/animal_lswallow_cycle_p_increase_std)

fprintf("\t\tModel Loaded Swallow Cycle Time: %.3f\n",model_data.lswallow_time)

metric_ind = metric_ind - 1;

%% III. b) Percent Protraction

model_lswallow_perc_prot = model_data.lswallow_prot_frac;
animal_lswallow_perc_prot_means = bootstrap_animal_data.lswallow_prot_frac_data;
animal_lswallow_time_std = bootstrap_animal_data.lswallow_prot_frac_std;
lswallow_perc_proc_frac_diff = (animal_lswallow_perc_prot_means - model_lswallow_perc_prot);
lswallow_perc_proc_90CI = quantile(lswallow_perc_proc_frac_diff,[0.05,0.95]);
lswallow_perc_proc_95CI = quantile(lswallow_perc_proc_frac_diff,[0.025,0.975]);

plot(d_lswallow*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot(mean(lswallow_perc_proc_frac_diff)/animal_lswallow_time_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(lswallow_perc_proc_95CI/animal_lswallow_time_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(lswallow_perc_proc_90CI/animal_lswallow_time_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

metric_ind = metric_ind - 1;

fprintf("\tPercent Protraction:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",mean(animal_lswallow_perc_prot_means),std(animal_lswallow_perc_prot_means),animal_lswallow_time_std)
fprintf("\t\tModel Value: %.3f\n",model_lswallow_perc_prot)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(mean(animal_lswallow_perc_prot_means) - model_lswallow_perc_prot)/mean(animal_lswallow_perc_prot_means), ...
                                                        (mean(animal_lswallow_perc_prot_means) - model_lswallow_perc_prot)/animal_lswallow_time_std)


%% IV. Rejection
fprintf("Rejection Metric Comparison:\n")

file = "Ye2006_RejectionTiming.csv";
data_reject_time = readmatrix("Animal Data\Datasets\"+file);
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



%% IV. a) Cycle Time 
model_reject_time = model_data.reject_time;
animal_reject_time_mean = t_end - t_i2_on;
animal_reject_time_std = sqrt(sig_end^2 + sig_i2_on^2);
animal_reject_90CI = (animal_reject_time_mean - model_reject_time) + tinv([0.05,0.95],5-1)*animal_reject_time_std/sqrt(5);
animal_reject_95CI = (animal_reject_time_mean - model_reject_time) + tinv([0.025,0.975],5-1)*animal_reject_time_std/sqrt(5);


plot(d_reject*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot((animal_reject_time_mean - model_reject_time)/animal_reject_time_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(animal_reject_95CI/animal_reject_time_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(animal_reject_90CI/animal_reject_time_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

metric_ind = metric_ind - 1;

fprintf("\tPercent Protraction:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",animal_reject_time_mean,animal_reject_time_std/sqrt(5),animal_reject_time_std)
fprintf("\t\tModel Value: %.3f\n",model_reject_time)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(animal_reject_time_mean - model_reject_time)/animal_reject_time_mean, ...
                                                        (animal_reject_time_mean - model_reject_time)/animal_reject_time_std)

%% IV. b) Percent Protraction
model_reject_perc_prot = model_data.reject_prot_frac;
animal_reject_perc_prot_mean = dur_i2 / animal_reject_time_mean;
animal_reject_perc_prot_std = animal_reject_perc_prot_mean*sqrt((sig_dur_i2/dur_i2)^2 + (animal_reject_time_std/animal_reject_time_mean)^2);
animal_reject_perc_prot_90CI = (animal_reject_perc_prot_mean - model_reject_perc_prot) + tinv([0.05,0.95],5-1)*animal_reject_perc_prot_std/sqrt(5);
animal_reject_perc_prot_95CI = (animal_reject_perc_prot_mean - model_reject_perc_prot) + tinv([0.025,0.975],5-1)*animal_reject_perc_prot_std/sqrt(5);


plot(d_reject*[-1,1],metric_ind*[1,1],"-",'Color',0.9*[1,1,1],"LineWidth",15)
plot((animal_reject_perc_prot_mean - model_reject_perc_prot)/animal_reject_perc_prot_std,metric_ind,'sk',"MarkerFaceColor","k","MarkerSize",8)
plot(animal_reject_perc_prot_95CI/animal_reject_perc_prot_std,[metric_ind,metric_ind],'-',"LineWidth",2,"Color",col_95)
plot(animal_reject_perc_prot_90CI/animal_reject_perc_prot_std,[metric_ind,metric_ind],'-k',"LineWidth",3)

metric_ind = metric_ind - 1;

fprintf("\tPercent Protraction:\n")
fprintf("\t\tAnimal Value: %.3f +/- %.3f [%.3f]\n",animal_reject_perc_prot_mean,animal_reject_perc_prot_std/sqrt(5),animal_reject_perc_prot_std)
fprintf("\t\tModel Value: %.3f\n",model_reject_perc_prot)
fprintf("\t\t%% Error: %.3f | Effect Size: %.3f\n",100*(animal_reject_perc_prot_mean - model_reject_perc_prot)/animal_reject_perc_prot_mean, ...
                                                        (animal_reject_perc_prot_mean - model_reject_perc_prot)/animal_reject_perc_prot_std)

%% Creating Labels
ylim([-0.5,7.5])
yticks(0:7)
yticklabels(labels)
xlabel("Standard Difference [Cohen's d]")
set(gca,"FontSize",15)