%% ~~~~~~~~~~~~~~~ Non-normalized Neustadter 2007 Comparison ~~~~~~~~~~~~~~~ 
%{
    This script generates Figure 7c, comparing the kinematics of bites and 
    unloaded swallows with the MRI data presented in Neustadter et al. 
    2007. Additional cosmetic corrections were then made in Inkscape, but 
    the data traces were unchanged. Data are realigned such that peak 
    protraction is aligned.
%}

t = linspace(0,1,100)'; % time vector to resample data

% finding the peak in the animal data
t_b_animal = t(animal_data.bite_trans_mean(t) == max(animal_data.bite_trans_mean(t)));          % point of peak protraction in biting
t_s_animal = t(animal_data.swallow_trans_mean(t) == max(animal_data.swallow_trans_mean(t)));    % point of peak protraction in swallowing
t_r_animal = t(animal_data.reject_trans_mean(t) == max(animal_data.reject_trans_mean(t)));      % point of peak protraction in rejection

%% Realinging Model Data with Determined Peak Protraction from Neustadter et al. 2007

if show_biting
    % Normalizing time to fraction of cycle time
    t_bite = (out_bite.tout(i_ss_bite)-t0b)/(out_bite.tout(i_ss_bite(end))-t0b);
    
    % converting model units to mm
    model_bite_trans = R_ref*10*( (out_bite.x_gh(i_ss_bite))*(0.5*(L0_I1v + L0_I1d)) - (0.5*(L0_I1v + L0_I1d)) );
    
    % converting radians to degree
    model_bite_angle = (180/pi)*(out_bite.theta_g_animal(i_ss_bite));
    
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
    model_bite_angle = interp1(t_bite,model_bite_angle,t);
    
    % reordering
    model_bite_trans = model_bite_trans(bite_reorder);
    model_bite_angle = model_bite_angle(bite_reorder);
end

if show_swallowing_unloaded
    % Normalizing time to fraction of cycle time
    t_uswallow = (out_uswallow.tout(i_ss_uswallow)-t0us)/(out_uswallow.tout(i_ss_uswallow(end))-t0us);

    % converting model units to mm
    model_swallow_trans = R_ref*10*( (out_uswallow.x_gh(i_ss_uswallow))*(0.5*(L0_I1v + L0_I1d))  - (0.5*(L0_I1v + L0_I1d)) );

    % converting radians to degree
    model_swallow_angle = (180/pi)*(out_uswallow.theta_g_animal(i_ss_uswallow));

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
    model_swallow_angle = interp1(t_uswallow,model_swallow_angle,t);
    
    % reordering
    model_swallow_trans = model_swallow_trans(swallow_reorder);
    model_swallow_angle = model_swallow_angle(swallow_reorder);
end

if show_rejection
    % Normalizing time to fraction of cycle time
    t_reject = (out_reject.tout(i_ss_reject)-t0r)/(out_reject.tout(i_ss_reject(end))-t0r);

    % converting model units to mm
    model_reject_trans = R_ref*10*( (out_reject.x_gh(i_ss_reject))*(0.5*(L0_I1v + L0_I1d))  - (0.5*(L0_I1v + L0_I1d)) );

    % converting radians to degree
    model_reject_angle = (180/pi)*(out_reject.theta_g_animal(i_ss_reject));

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
    model_reject_angle = interp1(t_reject,model_reject_angle,t);
    
    % reordering
    model_reject_trans = model_reject_trans(reject_reorder);
    model_reject_angle = model_reject_angle(reject_reorder);
end

%% Reading in Animal Data

% interpolation functions of digitized data
trans_bite_mean = @(t) animal_data.bite_trans_mean(t); 
trans_bite_std = @(t) animal_data.bite_trans_std(t);
trans_swallow_mean = @(t) animal_data.swallow_trans_mean(t);
trans_swallow_std = @(t) animal_data.swallow_trans_std(t); 
trans_reject_mean = @(t) animal_data.reject_trans_mean(t); 
trans_reject_std = @(t) 0*t; 

angle_bite_mean = @(t) animal_data.bite_angle_mean(t);
angle_bite_std = @(t) animal_data.bite_angle_std(t); 
angle_swallow_mean = @(t) animal_data.swallow_angle_mean(t); 
angle_swallow_std = @(t) animal_data.swallow_angle_std(t); 

%% Plotting
s = 3;
figure('Position',[0,0,2000,800],"Color","w");

err_range = [-80,80];
bound_range = 20*[-1,1];

if show_biting
% Biting Translation
    subplot(3,4,1); hold all
    
    t_ = [100*t; flipud(100*t)];
    err_ = [(trans_bite_mean(t) + trans_bite_std(t)); flipud((trans_bite_mean(t) - trans_bite_std(t)))];
    fill(t_,err_,animal_data_color,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");
    
    plot(100*t,(trans_bite_mean(t)),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Data');

    plot(100*t,model_bite_trans,'color',present_model_color,'LineWidth',linewidth,'DisplayName','Model')
    
    ylim([-11,3])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Biting","Translation [mm]"})
    legend("NumColumns",2)
    set(gca,"FontSize",15)

    % Biting Rotation
    subplot(3,4,2); hold all
    
    t_ = [100*t; flipud(100*t)];
    err_ = [(angle_bite_mean(t) + angle_bite_std(t)); flipud((angle_bite_mean(t) - angle_bite_std(t)))];
    fill(t_,err_,animal_data_color,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");
    
    plot(100*t,(angle_bite_mean(t) ),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Bite');
    plot(100*t,model_bite_angle,'color',present_model_color,'LineWidth',linewidth,'DisplayName','Model Bite')

            
    ylim([-10,115])
    ylabel({"Biting","Rotation [deg]"})
    xlabel('Normalized Time [%Cycle]')
    set(gca,"FontSize",15)

    subplot(3,4,3); hold all
    [max_b_t,max_b_t_ind] = max((trans_bite_mean(t)));
    [min_b_t,min_b_t_ind] = min((trans_bite_mean(t)));
    rom_b_t = max_b_t - min_b_t;
    std_rom_b_t = sqrt( (trans_bite_std(t(max_b_t_ind))).^2 + (trans_bite_std(t(min_b_t_ind))).^2 );
    

    err_b_t = ( trans_bite_mean(t) - model_bite_trans )/rom_b_t;
    sig_err_b_t = sqrt( ( trans_bite_std(t)./rom_b_t ).^2 + ( (err_b_t./rom_b_t).*std_rom_b_t ).^2 );

    t_ = [100*t; flipud(100*t)];
    err_ = 100*[(err_b_t + sig_err_b_t); flipud((err_b_t - sig_err_b_t))];
    fill(t_,err_,0.5*present_model_color,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");

    plot(100*t,100*err_b_t,'-','Color',present_model_color,'LineWidth',linewidth,'DisplayName','Animal Bite');
    
    plot([0:10:100],bound_range(1) + 0*[0:10:100],'--k')
    plot([0:10:100],bound_range(2) + 0*[0:10:100],'--k')
    
    
    ylim(err_range)
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Biting Translation","Error [%RoM]"})
    set(gca,"FontSize",15)

    subplot(3,4,4); hold all
    [max_b_r,max_b_r_ind] = max((angle_bite_mean(t)));
    [min_b_r,min_b_r_ind] = min((angle_bite_mean(t)));
    rom_b_r = max_b_r - min_b_r;
    std_rom_b_r = sqrt( (angle_bite_std(t(max_b_r_ind))).^2 + (angle_bite_std(t(min_b_r_ind))).^2 );
    

    err_b_r = ( angle_bite_mean(t) - model_bite_angle )/rom_b_r;
    sig_err_b_r = sqrt( ( angle_bite_std(t)./rom_b_r ).^2 + ( (err_b_r./rom_b_r).*std_rom_b_r ).^2 );

    t_ = [100*t; flipud(100*t)];
    err_ = 100*[(err_b_r + sig_err_b_r); flipud((err_b_r - sig_err_b_r))];
    fill(t_,err_,0.5*present_model_color,"EdgeColor","none","FaceAlpha",0.25,"HandleVisibility","off");    
    
    plot(100*t,100*err_b_r,'-','Color',present_model_color,'LineWidth',linewidth,'DisplayName','Animal Data');
    
    plot([0:10:100],bound_range(1) + 0*[0:10:100],'--k')
    plot([0:10:100],bound_range(2) + 0*[0:10:100],'--k')
    
    
    ylim(err_range)
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Biting Rotation","Error [%RoM]"})
    set(gca,"FontSize",15)

end

if show_swallowing_unloaded
    % Swallowing Translation
    subplot(3,4,5); hold all
    t_ = [100*t; flipud(100*t)];
    err_ = [(trans_swallow_mean(t) + trans_swallow_std(t)); flipud((trans_swallow_mean(t) - trans_swallow_std(t)))];
    fill(t_,err_,animal_data_color,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");

    plot(100*t,(trans_swallow_mean(t)),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Swallow');
    plot(100*t,model_swallow_trans,'color',present_model_color,'LineWidth',linewidth,'DisplayName','Model Swallow')
    
    ylim([-11,3])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Swallowing","Translation [mm]"})

    set(gca,"FontSize",15)
    
    % Swallowing Rotation
    subplot(3,4,6); hold all
    t_ = [100*t; flipud(100*t)];
    err_ = [(angle_swallow_mean(t) + angle_swallow_std(t)); flipud((angle_swallow_mean(t) - angle_swallow_std(t)))];
    fill(t_,err_,animal_data_color,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off");

    plot(100*t,(angle_swallow_mean(t) ),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Swallow');
    plot(100*t,model_swallow_angle,'color',present_model_color,'LineWidth',linewidth,'DisplayName','Model Swallow')
    
    ylim([-10,115])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Swallowing","Rotation [deg]"})
    set(gca,"FontSize",15)

    subplot(3,4,7); hold all
    [max_s_t,max_s_t_ind] = max((trans_swallow_mean(t)));
    [min_s_t,min_s_t_ind] = min((trans_swallow_mean(t)));
    rom_s_t = max_s_t - min_s_t;
    std_rom_s_t = sqrt( (trans_swallow_std(t(max_s_t_ind))).^2 + (trans_swallow_std(t(min_s_t_ind))).^2 );
    

    err_s_t = ( trans_swallow_mean(t) - model_swallow_trans )/rom_s_t;
    sig_err_s_t = sqrt( ( trans_swallow_std(t)./rom_s_t ).^2 + ( (err_s_t./rom_s_t).*std_rom_s_t ).^2 );
    
    t_ = [100*t; flipud(100*t)];
    err_ = 100*[(err_s_t + sig_err_s_t); flipud((err_s_t - sig_err_s_t))];
    fill(t_,err_,0.5*present_model_color,"EdgeColor","none","FaceAlpha",0.25,"HandleVisibility","off");
    
    plot(100*t,100*err_s_t,'-','Color',present_model_color,'LineWidth',linewidth,'DisplayName','Animal Bite');
    plot([0:10:100],bound_range(1) + 0*[0:10:100],'--k')
    plot([0:10:100],bound_range(2) + 0*[0:10:100],'--k')

    ylim(err_range)
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Swallowing Translation","Error [%RoM]"})
    set(gca,"FontSize",15)

    subplot(3,4,8); hold all
    [max_s_r,max_s_r_ind] = max((angle_swallow_mean(t)));
    [min_s_r,min_s_r_ind] = min((angle_swallow_mean(t)));
    rom_s_r = max_s_r - min_s_r;
    std_rom_s_r = sqrt( (angle_swallow_std(t(max_s_r_ind))).^2 + (angle_swallow_std(t(min_s_r_ind))).^2 );
    

    err_s_r = ( angle_swallow_mean(t) - model_swallow_angle )/rom_s_r;
    sig_err_s_r = sqrt( ( angle_swallow_std(t)./rom_s_r ).^2 + ( (err_s_r./rom_s_r).*std_rom_s_r ).^2 );

    t_ = [100*t; flipud(100*t)];
    err_ = 100*[(err_s_r + sig_err_s_r); flipud((err_s_r - sig_err_s_r))];
    fill(t_,err_,0.5*present_model_color,"EdgeColor","none","FaceAlpha",0.25,"HandleVisibility","off");
    
    plot(100*t,100*err_s_r,'-','Color',present_model_color,'LineWidth',linewidth,'DisplayName','Animal Bite');
    plot([0:10:100],bound_range(1) + 0*[0:10:100],'--k')
    plot([0:10:100],bound_range(2) + 0*[0:10:100],'--k')
    
    
    ylim(err_range)
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Swallowing Rotation","Error [%RoM]"})
    set(gca,"FontSize",15)
    
end


if show_rejection
    % Rejection Translation
    subplot(3,4,9); hold all
    plot(100*t,(trans_reject_mean(t)),'--','Color',animal_data_color,'LineWidth',linewidth,'DisplayName','Animal Rejection');
    plot(100*t,model_reject_trans,'color',present_model_color,'LineWidth',linewidth,'DisplayName','Model Rejection')
    
    t_ = [100*t; flipud(100*t)];

    ylim([-11,3])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Rejection","Translation [mm]"})
    
    set(gca,"FontSize",15)

    % Rejection Rotation
    subplot(3,4,10); hold all
    plot(100*t,model_reject_angle,'color',present_model_color,'LineWidth',linewidth,'DisplayName','Model Rejection')
    ylim([-10,115])

    % ylim([-10,110])
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Rejection","Rotation [deg]"})
    set(gca,"FontSize",15)

    subplot(3,4,11); hold all
    [max_r_t,max_r_t_ind] = max((trans_reject_mean(t)));
    [min_r_t,min_r_t_ind] = min((trans_reject_mean(t)));
    rom_r_t = max_r_t - min_r_t;
    std_rom_r_t = sqrt( (trans_reject_std(t(max_r_t_ind))).^2 + (trans_reject_std(t(min_r_t_ind))).^2 );
    

    err_r_t = ( trans_reject_mean(t) - model_reject_trans )/rom_r_t;
    sig_err_r_t = sqrt( ( trans_reject_std(t)./rom_r_t ).^2 + ( (err_r_t./rom_r_t).*std_rom_r_t ).^2 );

    plot(100*t,100*err_r_t,'-','Color',present_model_color,'LineWidth',linewidth,'DisplayName','Animal Bite');
    plot([0:10:100],bound_range(1) + 0*[0:10:100],'--k')
    plot([0:10:100],bound_range(2) + 0*[0:10:100],'--k')
    
    ylim(err_range)
    xlabel('Normalized Time [%Cycle]')
    ylabel({"Rejection Translation","Error [%RoM]"})
    set(gca,"FontSize",15)

end

%% Quantitative Comparisons

fprintf("Quantitative Non-normalized Comparisons for New Model\n")

[r,lags] = xcorr(model_bite_trans,trans_bite_mean(t),"normalized");
ind = find(lags==0);
R_bite_trans = r(ind);

RMSE_bite_trans = sqrt(mean((model_bite_trans - trans_bite_mean(t)).^2));

[r,lags] = xcorr(model_swallow_trans,trans_swallow_mean(t),"normalized");
ind = find(lags==0);
R_swallow_trans = r(ind);

RMSE_swallow_trans = sqrt(mean((model_swallow_trans - trans_swallow_mean(t)).^2));

[r,lags] = xcorr(model_reject_trans,trans_reject_mean(t),"normalized");
ind = find(lags==0);
R_reject_trans = r(ind);

RMSE_reject_trans = sqrt(mean((model_reject_trans - trans_reject_mean(t)).^2));

[r,lags] = xcorr(model_bite_angle,angle_bite_mean(t),"normalized");
ind = find(lags==0);
R_bite_angle = r(ind);

RMSE_bite_angle = sqrt(mean((model_bite_angle - angle_bite_mean(t)).^2));

[r,lags] = xcorr(model_swallow_angle,angle_swallow_mean(t),"normalized");
ind = find(lags==0);
R_swallow_angle = r(ind);

RMSE_swallow_angle = sqrt(mean((model_swallow_angle - angle_swallow_mean(t)).^2));

fprintf("\tCross-correlaion Translation\n");
fprintf("\t\tBiting: %.3f\n",R_bite_trans)
fprintf("\t\tSwallowing: %.3f\n",R_swallow_trans)
fprintf("\t\tRejection: %.3f\n",R_reject_trans)

fprintf("\tCross-correlaion Rotation\n");
fprintf("\t\tBiting: %.3f\n",R_bite_angle)
fprintf("\t\tSwallowing: %.3f\n",R_swallow_angle)

fprintf("\tRMSE Translation\n");
fprintf("\t\tBiting: %.3f mm\n",RMSE_bite_trans)
fprintf("\t\tSwallowing: %.3f mm\n",RMSE_swallow_trans)
fprintf("\t\tRejection: %.3f mm\n",RMSE_reject_trans)

fprintf("\tRMSE Rotation\n");
fprintf("\t\tBiting: %.3f deg\n",RMSE_bite_angle)
fprintf("\t\tSwallowing: %.3f deg\n",RMSE_swallow_angle)

%% Comparing Range of Motion

fprintf("\n\n~~~~Range of Motion Comparisons~~~~\n");
fprintf("\tTranslation:\n")
fprintf("\t\t\t\tModel: \tAnimal Data (mean +/- std)\n")
fprintf("\t\tBiting: %.3f \t %.3f+/-%.3f\n",model_data.bite_trans_rom, animal_data.bite_trans_rom_mean,animal_data.bite_trans_rom_std)
fprintf("\t\tSwallow: %.3f \t %.3f+/-%.3f\n",model_data.uswallow_trans_rom, animal_data.swallow_trans_rom_mean,animal_data.swallow_trans_rom_std)
fprintf("\t\tReject: %.3f \t %.3f+/-[ ]\n",model_data.reject_trans_rom, animal_data.reject_trans_rom_mean)

fprintf("\tRotation:\n")
fprintf("\t\t\t\tModel: \tAnimal Data (mean +/- std)\n")
fprintf("\t\tBiting: %.3f \t %.3f+/-%.3f\n",model_data.bite_angle_rom, animal_data.bite_angle_rom_mean,animal_data.bite_angle_rom_std)
fprintf("\t\tSwallow: %.3f \t %.3f+/-%.3f\n",model_data.uswallow_angle_rom, animal_data.swallow_angle_rom_mean,animal_data.swallow_angle_rom_std)

%% How much of the cycle duration was outside of the +/- 20% range
% not necessary for biting

% unloaded swallowing
ind_dur_over_s = find(abs(err_s_t) >= 0.2 );
t_dur_over_s = t(ind_dur_over_s);
dur_over_s = max(t_dur_over_s) - min(t_dur_over_s);
max_err_s = max(abs(err_s_t));

fprintf("\t\tIn swallowing, %.1f%% of the cycle was spent out of the +/-20%% error bound, with a max error %.1f%%.\n",dur_over_s*100,max_err_s*100)

% rejection
ind_dur_over_r = find(abs(err_r_t) >= 0.2 );
t_dur_over_r = t(ind_dur_over_r);
dur_over_r = max(t_dur_over_r) - min(t_dur_over_r);
max_err_r = max(abs(err_r_t));

fprintf("\t\tIn rejection, %.1f%% of the cycle was spent out of the +/-20%% error bound, with a max error %.1f%%.\n",dur_over_r*100,max_err_r*100)




