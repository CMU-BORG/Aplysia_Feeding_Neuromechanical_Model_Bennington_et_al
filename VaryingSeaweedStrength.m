%% Testing Response to Seaweed Strength
%{ 
    This scripts runs the simulated experiments presented in Figure 7a
    where the model animal feeds on seaweed of finite strength. The
    strength of the seaweed is swept for 0 to just stronger than the animal
    can break with 5 values in the sweep.
%}

% generating the strengths of seaweed
num_strengths = 5;
seaweed_strength_min = 0;
seaweed_strength_max = 0.25;
seaweed_strength_vec = linspace(seaweed_strength_min,seaweed_strength_max,num_strengths);

%       Setting the stimuli parameters
fixation_type = 0;
sens_mechanical_lips = 1;
sens_mechanical_grasper = 1;
sens_chemical_lips = 1;

CBI3_0 = 1;
B31B32_0 = 1;
B64_0 = 0;
B6B9B3_0 = 0;

t_switch_fixation = 4;   
t_switch_mech_lips = 1e6;
t_switch_chem_lips = 1e6;
t_switch_mech_grasp = 1e6;

tend = t_switch_fixation + 15;


figure("Position",[100,100,400,700]);
set(gcf,'Color','white')
for ind = 1:num_strengths
    % set the current seaweed_strength
    seaweed_strength = seaweed_strength_vec(ind);
    
    % run the simuation
    out = sim(simulinkFile); 
    % post process to calculate length of ingested seaweed
    CalculateLengthIngested;
    
    t = out.tout;
    % only get the points just before food becomes fixed
    pts = t>=t_switch_fixation-2;

    t = t(pts); t = t - t(1);
    F = out.F_fg(pts) + out.F_fh(pts);
    fixed = out.fixation(pts);  % indicator determining if the food is fixed and not broken
    t_fixed = t(fixed>0);       % find the region of time when the seaweed is fixed

    %% Plotting
    subplot(num_strengths,1,ind); hold all
    
    % shaded region showing when food is fixed
    fill([t_fixed(1),t_fixed(1),t_fixed(end),t_fixed(end)],[-0.2,2.2,2.2,-0.2],0.8*[1,1,1],"FaceAlpha",0.4,"EdgeColor","none")

    % plotting the normalized force
    plot(t,F / seaweed_strength_max,'k','LineWidth',2)

    % finding when the grasper is closed
    grasper_motion = out.x_gh(pts);
    grasper_pressure = out.P_I4(pts);
    idx = find(grasper_pressure >= 0.5);
    idy = find(grasper_pressure < 0.5);
    
    grasper_motion_pressure = 0*grasper_motion;

    grasper_motion_pressure(idx) = grasper_motion(idx);
    grasper_motion_pressure(idy)=NaN;

    % plotting the grasper motion (+1 allows the two traces to not overlap
    % on the same axis)
    plot(t,1+grasper_motion_pressure,'b','LineWidth',4)
    plot(t,1+grasper_motion,'b','LineWidth',2)
    hold off
    
    if i==1
        legend({'Normalized Force on Transducer','Normalized Grasper Motion'},'Orientation','horizontal','Position',[0.337905404191644,0.941469018401128,0.350331117341061,0.027292575735973],'Box','off','FontSize',12)
    end
    
    set(gca,'FontSize',15)
    grid on
    set(gca,'YGrid','off')
    ylim([-0.2 2.2])
    xlim([0,15])
    ylabel('Amplitude')
    set(gca,'ytick',[0 1]);
    set(gca,'YTickLabel',[]);
    
    i=i+1;
end
xlabel('Time (s)')
