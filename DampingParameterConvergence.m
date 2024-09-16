%% Dampering Parameter Convergence
%{
    This script tests the values what values are tau_I2 (the time scalar
    for the model) is required for behavior timing agreement as we decrease
    the value of the damping parameter (asymptotically approaching a fully
    quasistatic model).
%}

if recalc_data
    N_points = 20; % number of damping parameters to test
    
    inv_C = logspace(0,log10(333.33),N_points);
    t_cycle = zeros(N_points,1);
    tau_I2s = zeros(N_points,1);
        
    for kk=1:N_points
        % resetting the time scalar to 1
        time_scaler = 1;
    
        % updating the damping parameter
        C_xg = (1/inv_C(kk))*time_scaler*force_scaler;
    
        ParameterRecalculation; % update all other parameters that are affected by damping and time scalar
        MCC_0 = 1; % initializing metacerebral cell activity (1: model "hungry", 0: model "satiated") 
        tend = 40; % less time is required because it is known that steady state will occur in <5s from parameter tuning
    
        % running biting behavior
        behavior = "biting";
        
        % running simulation
        GetSteadyState; 
        
        % getting steady state cycle time
        t_cycle(kk) = out.tout(i_ss(end))-t0;
        
        % calculating I2 time constant from animal data and cycle time
        tau_I2s(kk) = animal_data.bite_time_mean / t_cycle(kk);
    
    end
    out_data = {inv_C,t_cycle,tau_I2s};
    save(damping_data_file,"out_data");

else
    % just load the precalculated data
    out_data = load(damping_data_file).out_data;
    inv_C = out_data{1};
    t_cycle = out_data{2};
    tau_I2s = out_data{3};
end

%% Plotting convergence against the inverse of damping coefficient
figure("Position",[100,100,600,400],"Color","w"); 
semilogx(1./inv_C,tau_I2s,'ok',"MarkerFaceColor","k")
hold all
xline(1/50,'--',"color",[1,1,1]*0.5,"LineWidth",1.5)
yline(1.1,'--b',"LineWidth",1.5)
yline(1.16,'--r',"LineWidth",1.5)

xlabel('Damping Coefficient [normalized]')
ylabel("Predicted \tau_{I2} [s]")

set(gca, "FontSize",15)