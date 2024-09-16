%% Plotting Unloaded Swallowing Length Measured from Model v Lum et al 2006

t = linspace(0,1,100)';

%% Animal Data
L_animal_mean = animal_data.swallow_length_mean(t);
L_animal_std = animal_data.swallow_length_std(t);

% finding minimum of swallowing trace in the animal for alignment with
% model data
[~,ia_min] = min(L_animal_mean);

trans_swallow_mean = @(t) animal_data.swallow_trans_mean(t); 
trans_swallow_std = @(t) animal_data.swallow_trans_std(t); 

%% Model Results 
mod_color = [0.8,0.2,0.6];

t_mod = out_uswallow.tout(i_ss_uswallow);
t_mod = (t_mod - t_mod(1))/(t_mod(end) - t_mod(1));
L_mod = out_uswallow.L_ingested(i_ss_uswallow);
L_mod = L_mod - L_mod(1);

L_mod = interp1(t_mod,L_mod,t);

[~,i_min_mod] = min(L_mod);

if(i_min_mod-ia_min > 0)
    L_mod = [L_mod(i_min_mod-ia_min:end); L_mod(end) + L_mod(1:i_min_mod-ia_min-1)]; % rearranging
else
    L_mod = [L_mod(ia_min-i_min_mod:end); L_mod(end) + L_mod(1:ia_min-i_min_mod-1)]; % rearranging
end

% %%%%% Code copied from Neustadter2007_Comparison_NonNormalized %%%%%
% Normalizing time to fraction of cycle time
t_uswallow = (out_uswallow.tout(i_ss_uswallow)-t0us)/(out_uswallow.tout(i_ss_uswallow(end))-t0us);

% converting model units to mm
model_swallow_trans = R_ref*10*( (out_uswallow.x_gh(i_ss_uswallow))*(0.5*(L0_I1v + L0_I1d))  - (0.5*(L0_I1v + L0_I1d)) );

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
model_swallow_trans_mod = model_swallow_trans(swallow_reorder);


%% Plotting Comparison

figure("Position",[100,100,500,300],"Color","w"); hold all

% model swallow

% animal data
t_err = [t; flipud(t)];
L_err = [L_animal_mean + L_animal_std; flipud(L_animal_mean - L_animal_std) ];
fill(100*t_err',L_err',[1,1,1]*0.5,"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off")
plot(100*t,L_animal_mean,'--',"Color",[0,0,0],"DisplayName","Animal Data")
plot(100*t,L_mod*(0.5*(L0_I1v + L0_I1d))*R_ref*10,"-","Color",mod_color,"LineWidth",1.5,"DisplayName",'Model')

xlabel('Normalized Time [%Cycle]')
ylabel({"Length of Ingested","Seaweed [mm]"})
legend("NumColumns",2)

set(gca,"FontSize",12)


%% Calculating RMSE and Cross-Correlation for the length of ingested seaweed

[r,lags] = xcorr(L_animal_mean,L_mod*(0.5*(L0_I1v + L0_I1d))*R_ref*10,"normalized");
ind = find(lags==0);
R_L_ingest = r(ind);

RMSE_L_ingest = sqrt(mean((L_mod*(0.5*(L0_I1v + L0_I1d))*R_ref*10 - L_animal_mean).^2));

fprintf("Ingested Seaweed Cross-Correlation and RMSE\n")
fprintf("\t\tR = %.3f\n",R_L_ingest)
fprintf("\t\tRMSE = %.3f\n",RMSE_L_ingest)