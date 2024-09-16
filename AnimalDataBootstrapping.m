%% Bootstrapping Population Averages

file = "Animal Data\AnimalTimingData.xlsx";
data = readtable(file);

behavior = "Loaded Swallowing";
behavior_means = data.CycleDurationMean(data.Behavior==behavior);
behavior_sem = data.CycleDurationSEM(data.Behavior==behavior);

if behavior == "Biting"
    mean_model = 4.36;
elseif behavior == "Unloaded Swallowing"
    mean_model = 4.83;
end

figure; hold all
bar(1:length(behavior_means),behavior_means)
errorbar(1:length(behavior_means),behavior_means,behavior_sem)

N_animals = length(behavior_means);
N_samps = 10000*N_animals;

mus = zeros(N_samps,1);

for i=1:N_samps
    animal_ind = randi(N_animals,N_animals,1); % which animal distributions to sample from

    means = zeros(N_animals,1);
    for j=1:N_animals
        means(j) = normrnd(behavior_means(animal_ind(j)),behavior_sem(animal_ind(j)),1);
    end
    mus(i) = mean(means);

end

figure; 
histogram(mus)
title(behavior)

fprintf("Population Average %s Duration: %.2f +/- %.2f\n",behavior,mean(mus),std(mus))

figure; 
histogram(mus - mean_model)
title(behavior+" Difference from Model")

