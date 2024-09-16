%% ~~~~~~~~~~~~~~~ All Behavior Comparison ~~~~~~~~~~~~~~~ 
%{
    This script generates Figure 5: Steady state simulation results using
    code modified from Webster-Wood et al. 2020. Each column shows a steady
    state cycle for a different behavioral condition. Additional cosmetic
    corrections were then made in Inkscape (including the addition of
    behavioral labels), but the data traces were unchanged.
%}

if_label = 1;
if_save = 1;

folder = output_folder; % location of the output data
files = {'BiteOutput.mat','UnloadedSwallowOutput.mat','LoadedSwallowOutput.mat','RejectOutput.mat'}; % output file names

outputfile = 'SuperFigure';
%% Figure Settings
labels = {'Biting','Unloaded Swallowing','Loaded Swallowing','Rejection'};
cols = length(labels);

% figure settings
figure('Position', [10 10 500+500*(cols/4) 1200]);
set(gcf,'Color','white')
xl=[0,1]; % show full time scale
ymin = 0;
ymax = 1;
shift = 0.0275;%0.0475;
top = 0.95;
i=0;
left = 0.25;
width = 0.7*(4/cols)/4.5;
hshift = 0.7/cols;
height = 0.015;
xticks([0:0.1:1])
subplot(24,cols,1)

data_labels = {"bite","uswallow","lswallow","reject"};

%% Loop through all of the behaviors to plot
for j=1:cols
    if_label = j==1; % y axis labels are only generated for the leftmost column
    i = 0;
    
    % reading in simulation data
    data_struct = load([folder,files{j}]).("out_"+data_labels{j});

    t = data_struct.tout;   % time vector of the simulation
    
    % finding steady state cycle
    dB31B32 = diff(data_struct.B31B32); % getting steady state cycles
    [~,starts] = findpeaks(dB31B32,"MinPeakDistance",1000);
    pts = starts(end-2):starts(end-1); % indexes of the steady state cycle
    
    % normalizing the time vector to fraction of cycle length
    t = t(pts);
    t = t - t(1);
    t = t/max(t);

    xlimits = [0, 1];

    % breaking out simulation structure into individual variables
    %   cerebral-buccal interconnects
    CBI2 = data_struct.CBI2(pts);
    CBI3 = data_struct.CBI3(pts);
    CBI4 = data_struct.CBI4(pts);

    % interneurons
    B64 = data_struct.B64(pts);
    B20 = data_struct.B20(pts);
    B40B30 = data_struct.B40B30(pts);
    B4B5 = data_struct.B4B5(pts);
    
    % motor neurons
    B31B32 = data_struct.B31B32(pts);
    B8 = data_struct.B8(pts);
    B38 = data_struct.B38(pts);
    B6B9B3 = data_struct.B6B9B3(pts);
    B7 = data_struct.B7(pts);
    B82 = data_struct.B82(pts);
    C1C6 = data_struct.C1(pts);
    C2 = data_struct.C2(pts);

    % pressure in the grasper
    P_G = data_struct.P_I4(pts);
    
    % sensory cues
    AROUSAL = data_struct.MCC(pts);
    MECH_G = data_struct.out_sens_mechanical_grasper(pts);
    CHEM_L = data_struct.out_sens_chemical_lips(pts);
    MECH_L = data_struct.out_sens_mechanical_lips(pts);
    fixed = data_struct.fixation(pts);
    
    % model observables
    x_gh = data_struct.x_gh(pts);
    theta_g = data_struct.theta_g_animal(pts);
    dL_s = data_struct.L_ingested(pts);
    dL_s = dL_s - dL_s(1);
    force = data_struct.F_fg(pts,1) + data_struct.F_fh(pts,1);


    %% External Stimuli
    subplot('position',[left+(j-1)*hshift top width height])
    hold all
    i=i+1;
    plot(t,MECH_G, 'Color', [56/255, 232/255, 123/255],'LineWidth',2) %mechanical in grasper
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([0 1]);
    grid on
    
    xlim(xl)
    if if_label
        ylabel('Mech. in Grasper')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])    
    end
    set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    i=i+1;
    plot(t,CHEM_L, 'Color', [70/255, 84/255, 218/255],'LineWidth',2) %chemical at lips
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([0 1])
    grid on
    xlim(xl)

    if if_label
        ylabel('Chem. at Lips')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
    set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    i=i+1;
    plot(t,MECH_L, 'Color', [47/255, 195/255, 241/255],'LineWidth',2) %mechanical at lips
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([0 1])
    grid on
    xlim(xl)

    if if_label
        ylabel('Mech. at Lips')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
        
    end
    set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    i=i+1;
    plot(t,fixed, 'Color', [150/255, 50/255, 150/255],'LineWidth',2) %is seaweed fixed
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([0 1])
    grid on
    xlim(xl)

    if if_label
        ylabel('Seaweed Fixed')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
        
    end
    set(gca,'XColor','none')

    %% Cerebral Neurons
    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,CBI2,'k','LineWidth',2) % CBI2
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)


    if if_label
        ylabel('CBI-2')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
    set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,CBI3,'k','LineWidth',2) % CBI3
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)


    if if_label
        ylabel('CBI-3')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,CBI4,'k','LineWidth',2) % CBI4
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('CBI-4')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    %% Interneurons
    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B64,'LineWidth',2, 'Color',[90/255, 131/255, 198/255]) % B64
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)


    if if_label
        ylabel('B64', 'Color',[90/255, 131/255, 198/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B20,'LineWidth',2, 'Color',[44/255, 166/255, 90/255]) % B20
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('B20', 'Color',[44/255, 166/255, 90/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])

    end
        set(gca,'XColor','none')


    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B40B30,'LineWidth',2, 'Color',[192/255, 92/255, 185/255]) % B40/B30
    i=i+1.5;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)
    if if_label
        ylabel('B40/B30', 'Color',[192/255, 92/255, 185/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height*2])
    hold all
    plot(t,B4B5,'LineWidth',2, 'Color', [51/255, 185/255, 135/255]) % B4/5
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1,2])
    set(gca,'YTickLabel',[]);
    ylim([ymin 2])
    xlim(xl)

    if if_label
        ylabel('B4/B5', 'Color', [51/255, 185/255, 135/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    %% motor neurons
    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B31B32,'LineWidth',2, 'Color', [220/255, 81/255, 81/255]) % I2 input
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('B31/B32','Color',[220/255, 81/255, 81/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B8,'LineWidth',2, 'Color', [213/255, 155/255, 196/255]) % B8a/b
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('B8a/b', 'Color', [213/255, 155/255, 196/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B38,'LineWidth',2, 'Color', [238/255, 191/255, 70/255]) % B38
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('B38', 'Color', [238/255, 191/255, 70/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B6B9B3,'LineWidth',2, 'Color', [90/255, 155/255, 197/255]) % B6/9/3
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('B3/B6/B9', 'Color', [90/255, 155/255, 197/255])
        
        set(get(gca,'ylabel'),'rotation',0) 
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B7,'LineWidth',2, 'Color', [56/255, 167/255, 182/255]) % B7
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)


    if if_label
        ylabel('B7', 'Color', [56/255, 167/255, 182/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')

    % Adding in the two new neurons
    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,B82,'LineWidth',2, 'Color', [56/255, 167/255, 215/255]) % B7
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('*B82', 'Color', [56/255, 167/255, 215/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')


    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,C1C6,'LineWidth',2, 'Color', [25/255, 200/255, 25/255]) % C1/C6
    i=i+1;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('*C1/C6', 'Color', [25/255, 200/255, 25/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')


    subplot('position',[left+(j-1)*hshift top-i*shift width height])
    hold all
    plot(t,C2,'LineWidth',2, 'Color', [75/255, 175/255, 25/255]) % C1/C6
    i=i+3;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([ymin ymax])
    xlim(xl)

    if if_label
        ylabel('*C2', 'Color', [75/255, 175/255, 25/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')



    %% Model Observables

    %Grasper Motion
    subplot('position',[left+(j-1)*hshift top-i*shift width height*3.5])
    hold all
    grasper_motion = x_gh;
    grasper_pressure = P_G;

    idx = grasper_pressure >=0.5; % determining if the grasper is closed
    idy = grasper_pressure <0.5;
    grasper_motion_pressure = zeros(length(idx),1);
    grasper_motion_pressure(idx) = grasper_motion(idx);
    grasper_motion_pressure(idy)=NaN;

    plot(t,grasper_motion,'b','LineWidth',2)
    plot(t,grasper_motion_pressure','b','LineWidth',4)
   
    i=i+2.5;
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    set(gca,'ytick',[0,1])
    set(gca,'YTickLabel',[]);
    ylim([0,1.2])
    xlim(xl)


    if if_label
        ylabel({'Grasper';'Motion'}, 'Color', [0/255, 0/255, 255/255])
        set(gca,'XColor','none')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')


    subplot('position',[left+(j-1)*hshift top-i*shift width height*3.5])
    hold all
    plot(t,theta_g,'k','LineWidth',2)
    yticks([min(theta_g), mean(theta_g), max(theta_g)])

    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    xlim(xl)
    set(gca,'YTickLabel',[]);


    if if_label
        ylabel('Angle', 'Color', [0/255, 0/255, 0/255])
        set(gca,'XColor','none')
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
        set(gca,'XColor','none')
    i = i+2.5;


    subplot('position',[left+(j-1)*hshift top-i*shift width height*3.5])
    hold all
    plot(t,P_G,'k','LineWidth',2)
    yticks([0 1])
    set(gca,'YTickLabel',[]);
    set(gca,'FontSize',16)
    set(gca,'xtick',[])
    xlim(xl)
    ylim([0,1])
    set(gca,'XColor','none')


    if if_label
        ylabel({"Grasper","Pressure"}, 'Color', [0/255, 0/255, 0/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end
    i = i+2.5;
    
    
    subplot('position',[left+(j-1)*hshift top-i*shift width height*3.5])
    hold all
    plot(t,force,'k','LineWidth',2)
    ylim([-.1,0.3])
    yticks([-1 0 1])
    yticklabels({'','0',''})
    set(gca,'FontSize',16)
    set(gca,'xtick',[])

    xlim(xl)
    set(gca,'XColor','none')

    if if_label
        
        ylabel('Force', 'Color', [0/255, 0/255, 0/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
    end

    i=i+2.5;

    subplot('position',[left+(j-1)*hshift top-i*shift width height*3.5])
    hold all
    plot(t,dL_s,'k','LineWidth',2)
    ylim([-1,1])
    yticks([-1 0 1])
    yticklabels({'','0',''})
    set(gca,'FontSize',16)
    set(gca,'xtick',[])

    xlim(xl)
    set(gca,'XColor','none')

    if if_label
        
        ylabel({'Ingested','Seaweed'}, 'Color', [0/255, 0/255, 0/255])
        hYLabel = get(gca,'YLabel');
        set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right','Position',get(hYLabel,'Position')-[0.05 0 0])
        set(gca,'XColor','none')
    end
    

    
end
if if_save
   savefig([folder,outputfile]); 
end
