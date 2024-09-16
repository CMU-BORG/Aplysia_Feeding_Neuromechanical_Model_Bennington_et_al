%% Computational Power Analysis
clear all; close all; clc

%% Simulating Power of Equivalence Test - fixed N, vary effect size
num_exp = 10000; % number of experiments to run
rng(0);

% setting the parameters of the distribution being compared to the model
% results. Because the test deals in normalized differences and effect
% size, the the values chosen here do not affect the power analysis
animal_mean = 1;
animal_std = 1;

% colors for plotting the results of different sample sizes
colors = {[162,250,163]/255, [79,117,155]/255, [87,31,78]/255, [0,0,0]};

figure("Position",[100,100,1000,400],"Color","w"); hold all
ylim([0,1.2])

% setting up the range of effect sizes to test
d_max = 2;
dd = 0.1:0.01*d_max:d_max;
xlim([0,d_max])

% setting the observed difference to 0
model_mean = animal_mean;

ind = 0;
% for all sample sizes available in the current animal dataset
for k=[5,7,12]
    ind = ind + 1;

    % preallocating memory for the values estimates of power
    mus = zeros(length(dd),1);
    sigs = mus;
    
    % setting the current sample size
    N = k;

    % for all equivalence intervals of interest
    for j=1:length(dd)
        % setting the current equivalence interval
        d = dd(j);
        
        % preallocating memory for the resampled power values
        mu_ = zeros(100,1);

        % prealloicating memory for the indicator variable (1: equivalence
        % test passed)
        success = zeros(num_exp,1);
        for i=1:num_exp
            % randomly sample animal data from the true animal distribution
            animal_sample = normrnd(animal_mean,animal_std,N,1);

            % perform the equivalence test based on TOST procedure
            t1 = ((mean(animal_sample)-model_mean)+d*std(animal_sample))/(std(animal_sample)/sqrt(N));
            t2 = (-(mean(animal_sample)-model_mean)+d*std(animal_sample))/(std(animal_sample)/sqrt(N));
            p1 = 1-tcdf(t1,N-1);
            p2 = 1-tcdf(t2,N-1);
            
            % record whether or not the test succeeded
            success(i) = (max([p1,p2]) < 0.05);
        end
        
%         % from the N_exp experiments that were run, resample from them to
%         % get the average power and the error on the average power
%         N_sub = 100; % number of sub experiments to resample
%         for i=1:N_sub
%             % randomly sample the experiments above with replacement
%             mu_(i) = mean(randsample(success,1000,true));
%         end
        % store the mean and standard deviation of the power
        mus(j) = mean(success);
%         sigs(j) = std(success);  
        
    end

    col =colors{ind};
%     err_dd = [dd,fliplr(dd)];
%     err_sig = [mus' + sigs', fliplr(mus' - sigs')];
%     fill(err_dd,err_sig,col,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");
    plot(dd,mus,'-',"Color",col,'LineWidth',2,"DisplayName","N="+num2str(N))
    
    % find the value of d that yields a mean power of 0.8
    err = @(d) (interp1(dd,mus,d) - 0.8)^2;
    effect_size = fminsearch(err,dd(find(mus>0.8,1,"first")));

    plot(effect_size,0.8,'xk',"MarkerSize",10,"LineWidth",2,"HandleVisibility","off")

    fprintf("For N = %i, use a d = %.3f\n",k,effect_size)

    drawnow()
end

yline(0.8,'--',"Color",0.6*[1,1,1],"HandleVisibility","off")

xlabel("Effect Size")
ylabel("Power")
% title("Achievable Power")
legend()
set(gca,"FontSize",15,"FontName","Arial")

saveas(gcf,"Results Figures\PowerAnalysis.svg","svg")
