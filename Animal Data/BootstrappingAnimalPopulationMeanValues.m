%% Bootstrapping Animal Population Mean Values
% clear all; close all; clc;

%% Set up
N_per_animal = 10000; % number of samples per animal to calculate
rng(0); % setting the random number generator seed for reproducibility

file = "Animal Data\AnimalTimingData.xlsx";
data = readtable(file);

animal_data = struct();

%% I. Biting
behavior = "Biting";
%% I. a) Cycle Duration
behavior_means = data.CycleDurationMean(data.Behavior==behavior);
behavior_sem = data.CycleDurationSEM(data.Behavior==behavior);

N_animals = length(behavior_means);
N_samps = N_per_animal*N_animals;

% model_mean = model_data.bite_time;

resamples = zeros(N_samps,N_animals);
mus = zeros(N_samps,1);
stds = zeros(N_samps,1);

for i=1:N_samps
    animal_ind = randi(N_animals,N_animals,1); % which animal distributions to sample from
    means = zeros(N_animals,1);
    for j=1:N_animals
        means(j) = normrnd(behavior_means(animal_ind(j)),behavior_sem(animal_ind(j)),1); % for each animal in the random sample, sample from their mean distribution
    end
    resamples(i,:) = means;
    mus(i) = mean(means);
    stds(i) = std(means);
end

figure("Position",[100,100,800,350],"Color","w");
subplot(1,3,1); hold all
bar(1:length(behavior_means),behavior_means,"FaceColor",[0.8,0.3,0.5],"DisplayName","Animal Data")
errorbar(1:length(behavior_means),behavior_means,behavior_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
ylabel("Biting Cycle Duration [s]")
xlabel("Animal ID")
legend()
ylim([0,10])

subplot(1,3,2)
histogram(mus,"FaceColor",[0.8,0.3,0.5],"EdgeAlpha",0.2)
xlabel("Mean [s]")
sgtitle("Biting Cycle Duration - Bootstrapped Population Mean Sampling Distribution")

subplot(1,3,3)
histogram(stds,"FaceColor",[0.8,0.3,0.5],"EdgeAlpha",0.2)
xlabel("STD [s]")


% Saving data and statistics
animal_data.bite_time_resamples = resamples;

animal_data.bite_time_data = mus;
animal_data.bite_time_mean = mean(mus);

animal_data.bite_time_std_data = stds;
animal_data.bite_time_std = mean(stds);

animal_data.bite_time_95CI = quantile(mus,[0.025,0.975]);
animal_data.bite_time_90CI = quantile(mus,[0.05,0.95]);

for i=1:3
    subplot(1,3,i); hold all
    set(gca,"FontName","Arial","FontSize",12)
end

saveas(gcf,"Results Figures\"+behavior+"_CycleDuration.svg","svg")

%% I. b) Protraction Fraction
duration_means = data.CycleDurationMean(data.Behavior==behavior);
duration_sem = data.CycleDurationSEM(data.Behavior==behavior);
protraction_means = data.ProtractionDurationMean(data.Behavior==behavior);
protraction_sem = data.ProtractionDurationSEM(data.Behavior==behavior);

resamples = zeros(N_samps,N_animals);
N_animals = length(duration_means);
N_samps = N_per_animal*N_animals;

mus = zeros(N_samps,1);
stds = zeros(N_samps,1);
for i=1:N_samps
    animal_ind = randi(N_animals,N_animals,1); % which animal distributions to sample from
    means = zeros(N_animals,1);
    for j=1:N_animals
        delta = 100;
        while delta>0
            cycle_dur = normrnd(duration_means(animal_ind(j)),duration_sem(animal_ind(j)),1); % for each animal in the random sample, sample from their mean distribution
            prot_dur = normrnd(protraction_means(animal_ind(j)),protraction_sem(animal_ind(j)),1);
            delta = prot_dur - cycle_dur; % removing samples that have longer protraction durations than cycle times (impossible)
        end
        means(j) = prot_dur / cycle_dur;
    end
    resamples(i,:) = means;
    mus(i) = mean(means);
    stds(i) = std(means);
end


figure("Position",[100,100,800,350],"Color","w");
subplot(1,3,1); hold all
bar(1:length(duration_means),duration_means,"FaceColor",[0.8,0.3,0.5],"DisplayName","Cycle Length")
errorbar(1:length(duration_means),duration_means,duration_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
bar(1:length(protraction_means),protraction_means,"FaceColor",0.6*[0.8,0.3,0.5],"DisplayName","Protraction Duration")
errorbar(1:length(protraction_means),protraction_means,protraction_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
legend("NumColumns",2,"Location","best")
ylabel("Biting Protraction And Cycle Duration [s]")
xlabel("Animal ID")
ylim([0,10])

subplot(1,3,2); hold all
histogram(mus,"FaceColor",[0.8,0.3,0.5],"EdgeAlpha",0.2)
xlabel("Mean [s/s]")

subplot(1,3,3)
histogram(stds,"FaceColor",[0.8,0.3,0.5],"EdgeAlpha",0.2)
xlabel("STD [s/s]")

sgtitle("Biting Protraction Fraction - Bootstrapped Distributions")
animal_data.bite_prot_frac_resamples = resamples;

animal_data.bite_prot_frac_data = mus;
animal_data.bite_prot_frac_mean = mean(mus);

animal_data.bite_prot_frac_std_data = stds;
animal_data.bite_prot_frac_std = mean(stds);

animal_data.bite_prot_frac_95CI = quantile(mus,[0.025,0.975]);
animal_data.bite_prot_frac_90CI = quantile(mus,[0.05,0.95]);

for i=1:3
    subplot(1,3,i); hold all
    set(gca,"FontName","Arial","FontSize",12)
end

saveas(gcf,"Results Figures\"+behavior+"_ProtractionFraction.svg","svg")

%% II. Unloaded Swallowing
behavior = "Unloaded Swallowing";
%% II. a) Cycle Duration
behavior_means = data.CycleDurationMean(data.Behavior==behavior);
behavior_sem = data.CycleDurationSEM(data.Behavior==behavior);

N_animals = length(behavior_means);
N_samps = N_per_animal*N_animals;


resamples = zeros(N_samps,N_animals);
mus = zeros(N_samps,1);
stds = zeros(N_samps,1);

for i=1:N_samps
    animal_ind = randi(N_animals,N_animals,1); % which animal distributions to sample from
    means = zeros(N_animals,1);
    for j=1:N_animals
        means(j) = normrnd(behavior_means(animal_ind(j)),behavior_sem(animal_ind(j)),1); % for each animal in the random sample, sample from their mean distribution
    end
    resamples(i,:) = means;
    mus(i) = mean(means);
    stds(i) = std(means);
end

animal_data.uswallow_time_resamples = means;

animal_data.uswallow_time_data = mus;
animal_data.uswallow_time_mean = mean(mus);

animal_data.uswallow_time_std_data = stds;
animal_data.uswallow_time_std = mean(stds);

animal_data.uswallow_time_95CI = quantile(mus,[0.025,0.975]);
animal_data.uswallow_time_90CI = quantile(mus,[0.05,0.95]);

figure("Position",[100,100,800,350],"Color","w");
subplot(1,3,1); hold all
bar(1:length(behavior_means),behavior_means,"FaceColor",[0.3,0.8,0.5])
errorbar(1:length(behavior_means),behavior_means,behavior_sem,"vertical","ok","MarkerFaceColor","k")
ylabel("Unloaded Swallowing Cycle Duration [s]")
xlabel("Animal ID")
ylim([0,10])

subplot(1,3,2)
histogram(mus,"FaceColor",[0.3,0.8,0.5],"EdgeAlpha",0.2)
xlabel("Mean [s]")
sgtitle("Mean Unloaded Swallowing Cycle Duration - Bootstrapped Distributions")

subplot(1,3,3)
histogram(stds,"FaceColor",[0.3,0.8,0.5],"EdgeAlpha",0.2)
xlabel("STD [s]")

for i=1:3
    subplot(1,3,i); hold all
    set(gca,"FontName","Arial","FontSize",12)
end

saveas(gcf,"Results Figures\"+behavior+"_CycleDuration.svg","svg")

%% II. b) Protraction Fraction
duration_means = data.CycleDurationMean(data.Behavior==behavior);
duration_sem = data.CycleDurationSEM(data.Behavior==behavior);
protraction_means = data.ProtractionDurationMean(data.Behavior==behavior);
protraction_sem = data.ProtractionDurationSEM(data.Behavior==behavior);


N_animals = length(duration_means);
N_samps = N_per_animal*N_animals;

resamples = zeros(N_samps,N_animals);
mus = zeros(N_samps,1);
stds = zeros(N_samps,1);

for i=1:N_samps
    animal_ind = randi(N_animals,N_animals,1); % which animal distributions to sample from
    means = zeros(N_animals,1);
    for j=1:N_animals
        delta = 100;
        while delta>0
            cycle_dur = normrnd(duration_means(animal_ind(j)),duration_sem(animal_ind(j)),1); % for each animal in the random sample, sample from their mean distribution
            prot_dur = normrnd(protraction_means(animal_ind(j)),protraction_sem(animal_ind(j)),1);
            delta = prot_dur - cycle_dur; % removing samples that have longer protraction durations than cycle times (impossible)
        end
        means(j) = prot_dur / cycle_dur;
    end
    resamples(i,:) = means;
    mus(i) = mean(means);
    stds(i) = std(means);
end

animal_data.uswallow_prot_frac_resamples = resamples;

animal_data.uswallow_prot_frac_data = mus;
animal_data.uswallow_prot_frac_mean = mean(mus);

animal_data.uswallow_prot_frac_std_data = stds;
animal_data.uswallow_prot_frac_std = mean(stds);

animal_data.uswallow_prot_frac_95CI = quantile(mus,[0.025,0.975]);
animal_data.uswallow_prot_frac_90CI = quantile(mus,[0.05,0.95]);

figure("Position",[100,100,800,350],"Color","w");
subplot(1,3,1); hold all
bar(1:length(duration_means),duration_means,"FaceColor",[0.3,0.8,0.5],"DisplayName","Cycle Length")
errorbar(1:length(duration_means),duration_means,duration_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
bar(1:length(protraction_means),protraction_means,"FaceColor",0.6*[0.3,0.8,0.5],"DisplayName","Protraction Duration")
errorbar(1:length(protraction_means),protraction_means,protraction_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
legend("NumColumns",2,"Location","best")
ylabel({"Unloaded Swallowing Protraction","And Cycle Duration [s]"})
xlabel("Animal ID")
ylim([0,10])

subplot(1,3,2)
histogram(mus,"FaceColor",[0.3,0.8,0.5],"EdgeAlpha",0.2)
xlabel("Mean [s/s]")

subplot(1,3,3)
histogram(stds,"FaceColor",[0.3,0.8,0.5],"EdgeAlpha",0.2)
xlabel("STD [s/s]")


sgtitle("Unloaded Swallowing Protraction Fraction - Bootstrapped Distributions")

for i=1:3
    subplot(1,3,i); hold all
    set(gca,"FontName","Arial","FontSize",12)
end

saveas(gcf,"Results Figures\"+behavior+"_ProtractionFraction.svg","svg")

%% III. Loaded Swallows
behavior1 = "Loaded Swallowing";
behavior2 = "Unloaded Swallowing";

%% III. a) Percent Increase in Cycle Duration
lswallow_time_mean = data.CycleDurationMean((data.Behavior==behavior1)&(data.DataSource=="Gill et al 2020"));
lswallow_time_sem = data.CycleDurationSEM((data.Behavior==behavior1)&(data.DataSource=="Gill et al 2020"));
uswallow_time_mean = data.CycleDurationMean((data.Behavior==behavior2)&(data.DataSource=="Gill et al 2020"));
uswallow_time_sem = data.CycleDurationSEM((data.Behavior==behavior2)&(data.DataSource=="Gill et al 2020"));

N_animals = length(lswallow_time_mean);
N_samps = N_per_animal*N_animals;

resamples = zeros(N_samps,N_animals);
mus = zeros(N_samps,1);
stds = zeros(N_samps,1);

for i=1:N_samps
    animal_ind = randi(N_animals,N_animals,1); % which animal distributions to sample from
    means = zeros(N_animals,1);
    for j=1:N_animals
        uswallow = normrnd(uswallow_time_mean(animal_ind(j)),uswallow_time_sem(animal_ind(j)),1); % for each animal in the random sample, sample from their mean distribution
        lswallow = normrnd(lswallow_time_mean(animal_ind(j)),lswallow_time_sem(animal_ind(j)),1);
        means(j) = (lswallow - uswallow)/uswallow;
    end
    resamples(i,:) = means;
    mus(i) = mean(means);
    stds(i) = std(means);
end

animal_data.lswallow_cycle_p_increase_resamples = resamples;

animal_data.lswallow_cycle_p_increase_data = mus;
animal_data.lswallow_cycle_p_increase_mean = mean(mus);

animal_data.lswallow_cycle_p_increase_std_data = stds;
animal_data.lswallow_cycle_p_increase_std = mean(stds);

animal_data.lswallow_cycle_p_increase_95CI = quantile(mus,[0.025,0.975]);
animal_data.lswallow_cycle_p_increase_90CI = quantile(mus,[0.05,0.95]);

figure("Position",[100,100,800,350],"Color","w");
subplot(1,3,1); hold all
bar(1:length(lswallow_time_mean),lswallow_time_mean,"FaceColor",[0.3,0.5,0.8],"DisplayName","Loaded Swallows")
errorbar(1:length(lswallow_time_mean),lswallow_time_mean,lswallow_time_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
bar(1:length(uswallow_time_mean),uswallow_time_mean,"FaceColor",0.6*[0.3,0.5,0.8],"DisplayName","Unloaded Swallows")
errorbar(1:length(uswallow_time_mean),uswallow_time_mean,uswallow_time_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
legend("NumColumns",2,"Location","best")
ylabel({"Unloaded and Loaded Swallowing","Cycle Duration [s]"})
xlabel("Animal ID")
ylim([0,10])

subplot(1,3,2)
histogram(100*mus,"FaceColor",[0.3,0.5,0.8],"EdgeAlpha",0.2)
xlabel("Mean [%]")

subplot(1,3,3)
histogram(100*stds,"FaceColor",[0.3,0.5,0.8],"EdgeAlpha",0.2)
xlabel("STD [%]")


sgtitle("Loaded Swallow Percent Cycle Increase [%] - Bootstrapped Distributions")

for i=1:3
    subplot(1,3,i); hold all
    set(gca,"FontName","Arial","FontSize",12)
end

saveas(gcf,"Results Figures\"+behavior1+"_PercentCycleIncrease.svg","svg")

%% III. b) Protraction Fraction
behavior = "Loaded Swallowing";
duration_means = data.CycleDurationMean(data.Behavior==behavior);
duration_sem = data.CycleDurationSEM(data.Behavior==behavior);
protraction_means = data.ProtractionDurationMean(data.Behavior==behavior);
protraction_sem = data.ProtractionDurationSEM(data.Behavior==behavior);


N_animals = length(duration_means);
N_samps = N_per_animal*N_animals;

resamples = zeros(N_samps,N_animals);
mus = zeros(N_samps,1);
stds = zeros(N_samps,1);

for i=1:N_samps
    animal_ind = randi(N_animals,N_animals,1); % which animal distributions to sample from
    means = zeros(N_animals,1);
    for j=1:N_animals
        delta = 100;
        while delta>0
            cycle_dur = normrnd(duration_means(animal_ind(j)),duration_sem(animal_ind(j)),1); % for each animal in the random sample, sample from their mean distribution
            prot_dur = normrnd(protraction_means(animal_ind(j)),protraction_sem(animal_ind(j)),1);
            delta = prot_dur - cycle_dur; % removing samples that have longer protraction durations than cycle times (impossible)
        end
        means(j) = prot_dur / cycle_dur;
    end
    resamples(i,:) = means;
    mus(i) = mean(means);
    stds(i) = std(means);
end

animal_data.lswallow_prot_frac_resamples = resamples;

animal_data.lswallow_prot_frac_data = mus;
animal_data.lswallow_prot_frac_mean = mean(mus);

animal_data.lswallow_prot_frac_std_data = stds;
animal_data.lswallow_prot_frac_std = mean(stds);

animal_data.lswallow_prot_frac_95CI = quantile(mus,[0.025,0.975]);
animal_data.lswallow_prot_frac_90CI = quantile(mus,[0.05,0.95]);

figure("Position",[100,100,800,350],"Color","w");
subplot(1,3,1); hold all
bar(1:length(duration_means),duration_means,"FaceColor",[0.3,0.5,0.8],"DisplayName","Cycle Length")
errorbar(1:length(duration_means),duration_means,duration_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
bar(1:length(protraction_means),protraction_means,"FaceColor",0.6*[0.3,0.5,0.8],"DisplayName","Protraction Duration")
errorbar(1:length(protraction_means),protraction_means,protraction_sem,"vertical","ok","MarkerFaceColor","k","HandleVisibility","off")
legend("NumColumns",2,"Location","best")
ylabel({"Loaded Swallowing Protraction","And Cycle Duration [s]"})
xlabel("Animal ID")
ylim([0,10])

subplot(1,3,2)
histogram(mus,"FaceColor",[0.3,0.5,0.8],"EdgeAlpha",0.2)
xlabel("Mean [s/s]")

subplot(1,3,3)
histogram(stds,"FaceColor",[0.3,0.5,0.8],"EdgeAlpha",0.2)
xlabel("STD [s/s]")

sgtitle("Loaded Swallowing Protraction Fraction - Bootstrapped Distributions")

for i=1:3
    subplot(1,3,i); hold all
    set(gca,"FontName","Arial","FontSize",12)
end

saveas(gcf,"Results Figures\"+behavior+"_ProtractionFraction.svg","svg")

%% IV. Saving Data

save("BootstrappedAnimalData.mat","animal_data")
