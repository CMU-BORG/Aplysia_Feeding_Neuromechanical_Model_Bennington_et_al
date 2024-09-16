%% ~~~~~~~~~~~~~~~ Normalized Neustadter 2007 Comparison ~~~~~~~~~~~~~~~ 
%{
    This script generates Figure 10, comparing the behavioral kinematics of 
    the animal (Neustadter et al. 2007) and the current model with the 
    model presented in Webster-Wood et al. 2020. Additional cosmetic 
    corrections were then made in Inkscape, but the data traces were 
    unchanged. Data are realigned such that peak protraction is aligned.
    Because there are no units to the kinematics in the Webster-Wood model,
    all traces are range-of-motion normalized for each behavior.
%}

t = linspace(0,1,80)'; % time vector to resample data
t_b_animal = t(animal_data.bite_trans_mean(t) == max(animal_data.bite_trans_mean(t)));          % point of peak protraction in biting
t_s_animal = t(animal_data.swallow_trans_mean(t) == max(animal_data.swallow_trans_mean(t)));    % point of peak protraction in swallowing
t_r_animal = t(animal_data.reject_trans_mean(t) == max(animal_data.reject_trans_mean(t)));      % point of peak protraction in rejection

% if all scale: calculate a single scaling for all of the behaviors, rather
% than scaling each behavior independently
all_scale = 0;

%% Realinging Model Data with Determined Peak Protraction from Neustadter et al. 2007

if show_biting

    % NEW MODEL DATA

    % Normalizing time to fraction of cycle time
    t_bite = (out_bite.tout(i_ss_bite)-t0b)/(out_bite.tout(i_ss_bite(end))-t0b);
    
    % converting model units to mm
    model_bite_trans = out_bite.x_gh(i_ss_bite);
    
    % finding time points that coincide with peak protraction in each behaviors
    t_b_model = t_bite(find(model_bite_trans==max(model_bite_trans),1));
    
    % Reordering sequence of data so that peak protraction in each behavior
    % will align with the time point of peak protraction observed in the
    % animal. As all analyzed cycles are steady state, this is equivalent to 
    % sampling the cycle with a time shift.
    if t_b_model > t_b_animal
        bite_reorder = [find(t>=( t_b_model - t_b_animal )); find(t<( t_b_model - t_b_animal ))];
    else
        bite_reorder = [find(t>=( 1 + t_b_model - t_b_animal )); find(t<(1 + t_b_model - t_b_animal ))];
    end
    
    % interpolating to resample at the same rate as interpolated Neustadter
    % data
    model_bite_trans = interp1(t_bite,model_bite_trans,t);
    
    % reordering
    model_bite_trans = model_bite_trans(bite_reorder);

    % getting range of motion for normalization
    min_bite_trans = min(model_bite_trans);
    max_bite_trans = max(model_bite_trans);

end

if show_swallowing_unloaded
    % Normalizing time to fraction of cycle time
    t_uswallow = (out_uswallow.tout(i_ss_uswallow)-t0us)/(out_uswallow.tout(i_ss_uswallow(end))-t0us);

    % converting model units to mm
    model_swallow_trans = out_uswallow.x_gh(i_ss_uswallow);

    % finding time points that coincide with peak protraction in each behaviors
    t_s_model = t_uswallow(find(model_swallow_trans==max(model_swallow_trans),1));
    
    % Reordering sequence of data so that peak protraction in each behavior
    % will align with the time point of peak protraction observed in the
    % animal. As all analyzed cycles are steady state, this is equivalent to 
    % sampling the cycle with a time shift.
    if t_s_model > t_s_animal
        swallow_reorder = [find(t>=( t_s_model - t_s_animal )); find(t<( t_s_model - t_s_animal ))];
    else
        swallow_reorder = [find(t>=( 1 + t_s_model - t_s_animal )); find(t<( 1 + t_s_model - t_s_animal ))];
    end
    
    % interpolating to resample at the same rate as interpolated Neustadter
    % data
    model_swallow_trans = interp1(t_uswallow,model_swallow_trans,t);
    
    % reordering
    model_swallow_trans = model_swallow_trans(swallow_reorder);

    % getting range of motion for normalization
    min_swallow_trans = min(model_swallow_trans);
    max_swallow_trans = max(model_swallow_trans);

end

if show_rejection
    % Normalizing time to fraction of cycle time
    t_reject = (out_reject.tout(i_ss_reject)-t0r)/(out_reject.tout(i_ss_reject(end))-t0r);

    % converting model units to mm
    model_reject_trans = out_reject.x_gh(i_ss_reject);

    % finding time points that coincide with peak protraction in each behaviors
    t_r_model = t_reject(find(model_reject_trans==max(model_reject_trans),1));
    
    % Reordering sequence of data so that peak protraction in each behavior
    % will align with the time point of peak protraction observed in the
    % animal. As all analyzed cycles are steady state, this is equivalent to 
    % sampling the cycle with a time shift.
    if t_r_model > t_r_animal
        reject_reorder = [find(t>=( t_r_model - t_r_animal )); find(t<( t_r_model - t_r_animal ))];
    else
        reject_reorder = [find(t>=( 1 + t_r_model - t_r_animal )); find(t<( 1 + t_r_model - t_r_animal ))];
    end
    
    % interpolating to resample at the same rate as interpolated Neustadter
    % data
    model_reject_trans = interp1(t_reject,model_reject_trans,t);
    
    % reordering
    model_reject_trans = model_reject_trans(reject_reorder);

    % getting range of motion for normalization
    min_reject_trans = min(model_reject_trans);
    max_reject_trans = max(model_reject_trans);
    
end

% Normalizing the data
if all_scale
    min_model = min([min_bite_trans, min_swallow_trans, min_reject_trans]);
    max_model = max([max_bite_trans, max_swallow_trans, max_reject_trans]);
    
    model_bite_trans = (model_bite_trans - min_model)/(max_model - min_model);
    model_swallow_trans = (model_swallow_trans - min_model)/(max_model - min_model);
    model_reject_trans = (model_reject_trans - min_model)/(max_model - min_model);
else
    model_bite_trans = (model_bite_trans - min_bite_trans)/(max_bite_trans - min_bite_trans);
    model_swallow_trans = (model_swallow_trans - min_swallow_trans)/(max_swallow_trans - min_swallow_trans);
    model_reject_trans = (model_reject_trans - min_reject_trans)/(max_reject_trans - min_reject_trans);

end

%% Boolean Model Normalization and Alignment

if show_biting

    % BOOLEAN MODEL DATA
    t_bool = 0:aplysia_bite.TimeStep:aplysia_bite.EndTime; % shared time vector for all behaviors

    % Normalizing time to fraction of cycle time
    t_bite_bool = (t_bool(i_ss_bite_bool) - t0b_bool)/( t_bool(i_ss_bite_bool(end)) - t0b_bool );
    
    % converting model units to mm
    model_bite_trans_bool = aplysia_bite.x_g(i_ss_bite_bool) - aplysia_bite.x_h(i_ss_bite_bool);
    
    % finding time points that coincide with peak protraction in each behaviors
    t_b_bool = t_bite_bool(find(model_bite_trans_bool==max(model_bite_trans_bool),1));
    
    % Reordering sequence of data so that peak protraction in each behavior
    % will align with the time point of peak protraction observed in the
    % animal. As all analyzed cycles are steady state, this is equivalent to 
    % sampling the cycle with a time shift.
    if t_b_bool > t_b_animal
        bite_reorder = [find(t>=( t_b_bool - t_b_animal )); find(t<( t_b_bool - t_b_animal ))];
    else
        bite_reorder = [find(t>=( 1 + t_b_bool - t_b_animal )); find(t<(1 + t_b_bool - t_b_animal ))];
    end
    
    % interpolating to resample at the same rate as interpolated Neustadter
    % data
    model_bite_trans_bool = interp1(t_bite_bool,model_bite_trans_bool,t);
    
    % reordering
    model_bite_trans_bool = model_bite_trans_bool(bite_reorder);

    % getting range of motion for normalization
    min_bite_trans_bool = min(model_bite_trans_bool);
    max_bite_trans_bool = max(model_bite_trans_bool);

end

if show_swallowing_unloaded
    % Normalizing time to fraction of cycle time
    t_uswallow_bool = (t_bool(i_ss_uswallow_bool) - t0us_bool)/( t_bool(i_ss_uswallow_bool(end)) - t0us_bool );

    % converting model units to mm
    model_swallow_trans_bool = aplysia_uswallow.x_g(i_ss_uswallow_bool) - aplysia_uswallow.x_h(i_ss_uswallow_bool);

    % finding time points that coincide with peak protraction in each behaviors
    t_s_bool = t_uswallow_bool(find(model_swallow_trans_bool==max(model_swallow_trans_bool),1));
    
    % Reordering sequence of data so that peak protraction in each behavior
    % will align with the time point of peak protraction observed in the
    % animal. As all analyzed cycles are steady state, this is equivalent to 
    % sampling the cycle with a time shift.
    if t_s_bool > t_s_animal
        swallow_reorder = [find(t>=( t_s_bool - t_s_animal )); find(t<( t_s_bool - t_s_animal ))];
    else
        swallow_reorder = [find(t>=( 1 + t_s_bool - t_s_animal )); find(t<( 1 + t_s_bool - t_s_animal ))];
    end
    
    % interpolating to resample at the same rate as interpolated Neustadter
    % data
    model_swallow_trans_bool = interp1(t_uswallow_bool,model_swallow_trans_bool,t);
    
    % reordering
    model_swallow_trans_bool = model_swallow_trans_bool(swallow_reorder);

    % getting range of motion for normalization
    min_swallow_trans_bool = min(model_swallow_trans_bool);
    max_swallow_trans_bool = max(model_swallow_trans_bool);

end

if show_rejection
    % Normalizing time to fraction of cycle time
    t_reject_bool = (t_bool(i_ss_reject_bool)-t0r_bool)/(t_bool(i_ss_reject_bool(end))-t0r_bool);

    % converting model units to mm
    model_reject_trans_bool = aplysia_reject.x_g(i_ss_reject_bool) - aplysia_reject.x_h(i_ss_reject_bool);

    % finding time points that coincide with peak protraction in each behaviors
    t_r_bool = t_reject_bool(find(model_reject_trans_bool==max(model_reject_trans_bool),1));
    
    % Reordering sequence of data so that peak protraction in each behavior
    % will align with the time point of peak protraction observed in the
    % animal. As all analyzed cycles are steady state, this is equivalent to 
    % sampling the cycle with a time shift.
    if t_r_bool > t_r_animal
        reject_reorder = [find(t>=( t_r_bool - t_r_animal )); find(t<( t_r_bool - t_r_animal ))];
    else
        reject_reorder = [find(t>=( 1 + t_r_bool - t_r_animal )); find(t<( 1 + t_r_bool - t_r_animal ))];
    end
    
    % interpolating to resample at the same rate as interpolated Neustadter
    % data
    model_reject_trans_bool = interp1(t_reject_bool,model_reject_trans_bool,t);
    
    % reordering
    model_reject_trans_bool = model_reject_trans_bool(reject_reorder);

    % getting range of motion for normalization
    min_reject_trans_bool = min(model_reject_trans_bool);
    max_reject_trans_bool = max(model_reject_trans_bool);
    
end

% Normalizing the data
if all_scale
    min_bool = min([min_bite_trans_bool, min_swallow_trans_bool, min_reject_trans_bool]);
    max_bool = max([max_bite_trans_bool, max_swallow_trans_bool, max_reject_trans_bool]);
    
    model_bite_trans_bool = (model_bite_trans_bool - min_bool)/(max_bool - min_bool);
    model_swallow_trans_bool = (model_swallow_trans_bool - min_bool)/(max_bool - min_bool);
    model_reject_trans_bool = (model_reject_trans_bool - min_bool)/(max_bool - min_bool);
else
    model_bite_trans_bool = (model_bite_trans_bool - min_bite_trans_bool)/(max_bite_trans_bool - min_bite_trans_bool);
    model_swallow_trans_bool = (model_swallow_trans_bool - min_swallow_trans_bool)/(max_swallow_trans_bool - min_swallow_trans_bool);
    model_reject_trans_bool = (model_reject_trans_bool - min_reject_trans_bool)/(max_reject_trans_bool - min_reject_trans_bool);

end

%% Reading in Animal Data

% interpolation functions of digitized data
trans_bite_mean = animal_data.bite_trans_mean(t);
trans_bite_std = animal_data.bite_trans_std(t); 
trans_swallow_mean = animal_data.swallow_trans_mean(t);
trans_swallow_std = animal_data.swallow_trans_std(t); 
trans_reject_mean = animal_data.reject_trans_mean(t); 
trans_reject_std =  0*t; 


max_bite_trans_animal = max(trans_bite_mean);
min_bite_trans_animal = min(trans_bite_mean);

max_swallow_trans_animal = max(trans_swallow_mean);
min_swallow_trans_animal = min(trans_swallow_mean);

max_reject_trans_animal = max(trans_reject_mean);
min_reject_trans_animal = min(trans_reject_mean);

% Normalizing the animal data
if all_scale
    min_animal = min([ min_bite_trans_animal, min_swallow_trans_animal, min_reject_trans_animal ]);
    max_animal = max([ max_bite_trans_animal, max_swallow_trans_animal, max_reject_trans_animal ]);
    
    trans_bite_mean = (trans_bite_mean - min_animal)/(max_animal - min_animal);
    trans_swallow_mean = (trans_swallow_mean - min_animal)/(max_animal - min_animal);
    trans_reject_mean = (trans_reject_mean - min_animal)/(max_animal - min_animal);
    
    trans_bite_std = (trans_bite_std )/(max_animal - min_animal);
    trans_swallow_std = (trans_swallow_std)/(max_animal - min_animal);
else
    trans_bite_mean = (trans_bite_mean - min_bite_trans_animal)/(max_bite_trans_animal - min_bite_trans_animal);
    trans_swallow_mean = (trans_swallow_mean - min_swallow_trans_animal)/(max_swallow_trans_animal - min_swallow_trans_animal);
    trans_reject_mean = (trans_reject_mean - min_reject_trans_animal)/(max_reject_trans_animal - min_reject_trans_animal);
    
    trans_bite_std = (trans_bite_std )/(max_bite_trans_animal - min_bite_trans_animal);
    trans_swallow_std = (trans_swallow_std)/(max_swallow_trans_animal - min_swallow_trans_animal);
end


%% Plotting
s = 3;
figure('Position',[100,100,535,800],"Color","w");

err_range = [-80,80];

if show_biting
% Biting Translation
    subplot(3,1,1); hold all

    t_ = [100*t; flipud(100*t)];
    err_ = [(trans_bite_mean + trans_bite_std); flipud((trans_bite_mean - trans_bite_std))];
    fill(t_,100*err_,2*animal_data_color,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");

    plot(100*t,100*(trans_bite_mean),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Data');
    plot(100*t,100*model_bite_trans,'-','color',present_model_color,'LineWidth',linewidth,'DisplayName','Present Model')
    plot(100*t,100*model_bite_trans_bool,'-.','color',boolean_model_color,'LineWidth',linewidth,'DisplayName','Previous Model')
    
    ylim(100*[-0.1,1.1])
    xlim(100*[0,1])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Normalized Biting","Translation [%RoM]"})
    legend("NumColumns",3)
    set(gca,"FontSize",15)

end

if show_swallowing_unloaded
    % Swallowing Translation
    subplot(3,1,2); hold all
    t_ = [100*t; flipud(100*t)];
    err_ = 100*[(trans_swallow_mean + trans_swallow_std); flipud((trans_swallow_mean - trans_swallow_std))];
    fill(t_,err_,animal_data_color,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");

    plot(100*t,100*(trans_swallow_mean),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Swallow');
    plot(100*t,100*model_swallow_trans,'-','color',present_model_color,'LineWidth',linewidth,'DisplayName','Model Swallow')
    plot(100*t,100*model_swallow_trans_bool,'-.','color',boolean_model_color,'LineWidth',linewidth,'DisplayName','Original Model')

    ylim(100*[-0.1,1.1])
    xlim(100*[0,1])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Normalized Swallowing","Translation [%RoM]"})
    set(gca,"FontSize",15)
    
end


if show_rejection
    % Rejection Translation
    subplot(3,1,3); hold all

    plot(100*t,100*(trans_reject_mean),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Rejection');
    plot(100*t,100*model_reject_trans,'-','color',present_model_color,'LineWidth',linewidth,'DisplayName','Model Rejection')
    plot(100*t,100*model_reject_trans_bool,'-.','color',boolean_model_color,'LineWidth',linewidth,'DisplayName','Original Model')

    ylim(100*[-0.1,1.1])
    xlim(100*[0,1])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Normalized Rejection","Translation [%RoM]"})
    set(gca,"FontSize",15)
    
end

%% Quantitative Comparison
fprintf("Quantitative Comparisons between Original Boolean Model and New Model\n")

[r,lags] = xcorr(model_bite_trans,trans_bite_mean,"coef");
ind = find(lags==0);
R_bite_new = r(ind);

RMSE_bite_new = 100*sqrt(mean((model_bite_trans - trans_bite_mean).^2));

[r,lags] = xcorr(model_swallow_trans,trans_swallow_mean,"coef");
ind = find(lags==0);
R_swallow_new = r(ind);

RMSE_swallow_new = 100*sqrt(mean((model_swallow_trans - trans_swallow_mean).^2));

[r,lags] = xcorr(model_reject_trans,trans_reject_mean,"coef");
ind = find(lags==0);
R_reject_new = r(ind);

RMSE_reject_new = 100*sqrt(mean((model_reject_trans - trans_reject_mean).^2));

[r,lags] = xcorr(model_bite_trans_bool,trans_bite_mean,"coef");
ind = find(lags==0);
R_bite_bool = r(ind);

RMSE_bite_bool = 100*sqrt(mean((model_bite_trans_bool - trans_bite_mean).^2));

[r,lags] = xcorr(model_swallow_trans_bool,trans_swallow_mean,"coef");
ind = find(lags==0);
R_swallow_bool = r(ind);

RMSE_swallow_bool = 100*sqrt(mean((model_swallow_trans_bool - trans_swallow_mean).^2));

[r,lags] = xcorr(model_reject_trans_bool,trans_reject_mean,"coef");
ind = find(lags==0);
R_reject_bool = r(ind);

RMSE_reject_bool = 100*sqrt(mean((model_reject_trans_bool - trans_reject_mean).^2));

fprintf("\tCross-Correlation Values for Normalized Curves\n");
fprintf("\t\tBiting:\n\t\t\tNew Model: %.3f\n\t\t\tOriginal Model: %.3f\n",R_bite_new,R_bite_bool)
fprintf("\t\tSwallowing:\n\t\t\tNew Model: %.3f\n\t\t\tOriginal Model: %.3f\n",R_swallow_new,R_swallow_bool)
fprintf("\t\tRejection:\n\t\t\tNew Model: %.3f\n\t\t\tOriginal Model: %.3f\n",R_reject_new,R_reject_bool)

fprintf("\tRMSE Values for Normalized Curves\n");
fprintf("\t\tBiting:\n\t\t\tNew Model: %.3f\n\t\t\tOriginal Model: %.3f\n",RMSE_bite_new,RMSE_bite_bool)
fprintf("\t\tSwallowing:\n\t\t\tNew Model: %.3f\n\t\t\tOriginal Model: %.3f\n",RMSE_swallow_new,RMSE_swallow_bool)
fprintf("\t\tRejection:\n\t\t\tNew Model: %.3f\n\t\t\tOriginal Model: %.3f\n",RMSE_reject_new,RMSE_reject_bool)

% Timing Ratio Comparison
bite_swallow_animal = animal_data.unloaded_swallow_time_mean ./ animal_data.bite_time_mean;
dbite_swallow_animal = bite_swallow_animal*sqrt( (animal_data.unloaded_swallow_time_std / animal_data.unloaded_swallow_time_mean)^2 + (animal_data.bite_time_std / animal_data.bite_time_mean)^2 );
bite_swallow_model = (out_uswallow.tout(i_ss_uswallow(end))-t0us)/(out_bite.tout(i_ss_bite(end))-t0b);
bite_swallow_bool = (t_bool(i_ss_uswallow_bool(end)) - t0us_bool)/(t_bool(i_ss_bite_bool(end)) - t0b_bool);

bite_reject_animal = animal_data.reject_time_mean ./ animal_data.bite_time_mean;
dbite_reject_animal = bite_reject_animal*sqrt( (animal_data.reject_time_std / animal_data.reject_time_mean)^2 + (animal_data.bite_time_std / animal_data.bite_time_mean)^2 );
bite_reject_model = (out_reject.tout(i_ss_reject(end))-t0r)/(out_bite.tout(i_ss_bite(end))-t0b);
bite_reject_bool = (t_bool(i_ss_reject_bool(end)) - t0r_bool)/(t_bool(i_ss_bite_bool(end)) - t0b_bool);


fprintf("\tRatio of Behavior Cycle Duration\n");
fprintf("\t\tSwallowing : Biting\n")
fprintf("\t\t\tAnimal: %.3f+/-%.3f\n\t\t\tNew Model: %.3f (%.1f %%, %.2f STD)\n\t\t\tOriginal Model: %.3f (%.1f, %.2f STD)\n", ...
                        bite_swallow_animal,dbite_swallow_animal, ...
                        bite_swallow_model,100*(bite_swallow_model - bite_swallow_animal)/bite_swallow_animal,(bite_swallow_model - bite_swallow_animal)/dbite_swallow_animal, ...
                        bite_swallow_bool,100*(bite_swallow_bool - bite_swallow_animal)/bite_swallow_animal,(bite_swallow_bool - bite_swallow_animal)/dbite_swallow_animal);

fprintf("\t\tRejection : Biting\n")
fprintf("\t\t\tAnimal: %.3f+/-%.3f\n\t\t\tNew Model: %.3f (%.1f %%, %.2f STD)\n\t\t\tOriginal Model: %.3f (%.1f, %.2f STD)\n", ...
                        bite_reject_animal,dbite_reject_animal, ...
                        bite_reject_model,100*(bite_reject_model - bite_reject_animal)/bite_reject_animal,(bite_reject_model - bite_reject_animal)/dbite_reject_animal, ...
                        bite_reject_bool,100*(bite_reject_bool - bite_reject_animal)/bite_reject_animal,(bite_reject_bool - bite_reject_animal)/dbite_reject_animal);



