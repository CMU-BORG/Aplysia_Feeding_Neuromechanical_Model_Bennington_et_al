%% ~~~~~~~~~~~~~~~ Gill 2020 Comparison ~~~~~~~~~~~~~~~ 
%{
    This script generates Figure 7b, which compares simulated loaded 
    swallows with data presented in Gill and Chiel 2020. Additional 
    cosmetic corrections were then made in Inkscape, but the data traces 
    were unchanged.
%}


% Figure settings
rot = 0;
cols = 2;
LW = 1.5;
figure('Position',[60,60,400*cols,500],"Color",[1,1,1])

for i=1:2
    
    % Reading in digitized animal data
    if i==2
        animal_data_file = "Animal Data\Gill2020LoadedSwallow.mat";
    else
        animal_data_file = "Animal Data\Gill2020UnloadedSwallow.mat";
    end
    animal_data_loaded_swallow = load(animal_data_file).data_interp;
    
    
    
    t_animal = 100*linspace(0,1,1000);
    I2_start_animal = find(animal_data_loaded_swallow{2,1}{1,1} > 1,1,"first"); % find when I2 activity starts in the animal
    t_start_animal = t_animal(I2_start_animal);

    % Determining the start of steady state cycles with B38, as plotted in Gill
    % 2020
    
    if i==2
        dB38 = diff(out_lswallow.B38);
        dB31B32 = diff(out_lswallow.B31B32); % getting steady state cycles
        starts = find(dB31B32==1);
       
        i_ss1 = starts(end-2):starts(end-1); % indices of steady state cycles
        i_ss2 = starts(end-1):starts(end);
        % these are the indicies of two consecutive steady state cycles

        dB381 = diff(out_lswallow.B38(i_ss1));
        ind_start = find(dB381==1,1,'last'); ind_start = i_ss1(ind_start);
        dB382 = diff(out_lswallow.B38(i_ss2)); 
        ind_stop = find(dB382==1,1,'last'); ind_stop  = i_ss2(ind_stop);
        
        i_ss_swallow = [i_ss1(1:end-1),i_ss2]; 
        
    else
        dB38 = diff(out_uswallow.B38);
        dB31B32 = diff(out_uswallow.B31B32); % getting steady state cycles
        starts = find(dB31B32==1);
       
        i_ss1 = starts(end-3):starts(end-2); % indices of steady state cycles
        i_ss2 = starts(end-2):starts(end-1);
        % these are the indicies of two consecutive steady state cycles

        dB381 = diff(out_uswallow.B38(i_ss1));
        ind_start = find(dB381==1,1,'last'); ind_start = i_ss1(ind_start);
        dB382 = diff(out_uswallow.B38(i_ss2)); 
        ind_stop = find(dB382==1,1,'last'); ind_stop  = i_ss2(ind_stop);
        
        i_ss_swallow = [i_ss1(1:end-1),i_ss2]; 
        
    end
    
    
    % Normalizing time data to %Cycle Time
    if i==2
        T = out_lswallow.tout(i_ss_swallow);
        T0 = out_lswallow.tout(i_ss2(1)); % realigning so start of I2 activity align
    else
        T = out_uswallow.tout(i_ss_swallow);
        T0 = out_uswallow.tout(i_ss2(1));
    end
    T = T - T0;
    T = 100*T/T(end) + t_start_animal;

    
    %% B38 Neural Activity
    subplot(6,cols,i); hold all
    T_sample = linspace(0,100,1000);
    % Model Data
    if i==2
        plot(T,out_lswallow.B38(i_ss_swallow),'Color',[240,194,75]/255,'LineWidth',LW)
        B38_model = interp1(T,out_lswallow.B38(i_ss_swallow),t_animal,"linear",0);
    else
        plot(T,out_uswallow.B38(i_ss_swallow),'Color',[240,194,75]/255,'LineWidth',LW)
        B38_model = interp1(T,out_uswallow.B38(i_ss_swallow),t_animal,"linear",0);
    end
    ylim([-0.1,1.5])
    if i==1
        ylabel("B38","Rotation",rot,"HorizontalAlignment","right","VerticalAlignment","middle","Color",'k')
    end
    yticks([0,1])
    set(gca,'YColor','k')
    
    % Animal Data
    yyaxis right
    hold all

    B38_animal = animal_data_loaded_swallow{1,1}{1,1};
    B38_animal_upper = animal_data_loaded_swallow{1,1}{2,1};
    B38_animal_lower = animal_data_loaded_swallow{1,1}{3,1};

    fill([t_animal,fliplr(t_animal)],[B38_animal_lower,fliplr(B38_animal_upper)],0.8*[240,194,75]/255,"EdgeColor","none","FaceAlpha",0.4)
    plot(t_animal,B38_animal,'-','Color',0.6*[240,194,75]/255,'LineWidth',LW)

    [r,lags] = xcorr(B38_animal,B38_model,"normalized");
    ind = find(lags==0);
    R_B38 = r(ind);

    ylim([-0.5,20.5])
    yticks([])

    text(100,10.25,"R = "+num2str(R_B38,2))
    
    xlim([0,100])
    xticks([])
    
    %% I2 Activity
    subplot(6,cols,i+2)
    
    % Model Data
    if i==2
        plot(T,out_lswallow.A_I2(i_ss_swallow),'color',[220,103,102]/255,'LineWidth',LW)
        I2_model = interp1(T,out_lswallow.A_I2(i_ss_swallow),t_animal,"linear",0);
    else
        plot(T,out_uswallow.A_I2(i_ss_swallow),'color',[220,103,102]/255,'LineWidth',LW)
        I2_model = interp1(T,out_uswallow.A_I2(i_ss_swallow),t_animal,"linear",0);
    end
    if i==1
        ylabel("I2 Activity","Rotation",rot,"HorizontalAlignment","right","VerticalAlignment","middle","Color",'k')
    end
    ylim([-0.1,1.1])
    yticks([0,1])
    set(gca,'YColor','k')
    
    % Animal Data
    yyaxis right
    hold all
    I2_animal = animal_data_loaded_swallow{2,1}{1,1};
    I2_animal_upper = animal_data_loaded_swallow{2,1}{2,1};
    I2_animal_lower = animal_data_loaded_swallow{2,1}{3,1};

    fill([t_animal,fliplr(t_animal)],[I2_animal_lower,fliplr(I2_animal_upper)],0.8*[220,103,102]/255,"EdgeColor","none","FaceAlpha",0.4)
    plot(t_animal,I2_animal,'-','Color',0.6*[220,103,102]/255,'LineWidth',LW)

    [r,lags] = xcorr(I2_animal,I2_model,"normalized");
    ind = find(lags==0);
    R_I2 = r(ind);

    yticks([])

    text(100,10,"R = "+num2str(R_I2,2))

    ylim([-2,22])
    yticks([])
    
    xlim([0,100])
    xticks([])
    
    %% B8 Activity
    subplot(6,cols,i+4)
    
    % Model Data
    if i==2
        plot(T,out_lswallow.B8(i_ss_swallow),'Color',[218,141,192]/255,'LineWidth',LW)
        B8_model = interp1(T,out_lswallow.B8(i_ss_swallow),t_animal,"linear",0);
    else
        plot(T,out_uswallow.B8(i_ss_swallow),'Color',[218,141,192]/255,'LineWidth',LW)
        B8_model = interp1(T,out_uswallow.B8(i_ss_swallow),t_animal,"linear",0);
    end
    if i==1
        ylabel("B8a/b","Rotation",rot,"HorizontalAlignment","right","VerticalAlignment","middle","Color",'k')
    end
    ylim([-0.1,1.1])
    yticks([0,1])
    set(gca,'YColor','k')
    
    % Animal Data
    yyaxis right
    hold all
    B8_animal = animal_data_loaded_swallow{3,1}{1,1};
    B8_animal_upper = animal_data_loaded_swallow{3,1}{2,1};
    B8_animal_lower = animal_data_loaded_swallow{3,1}{3,1};

    fill([t_animal,fliplr(t_animal)],[B8_animal_lower,fliplr(B8_animal_upper)],0.8*[218,141,192]/255,"EdgeColor","none","FaceAlpha",0.4)
    plot(t_animal,B8_animal,'-','Color',0.6*[218,141,192]/255,'LineWidth',LW)

    [r,lags] = xcorr(B8_animal,B8_model,"normalized");
    ind = find(lags==0);
    R_B8 = r(ind);

    yticks([])

    text(100,22,"R = "+num2str(R_B8,2))

    ylim([-0.1,1.1]*40)
    yticks([])
    
    xlim([0,100])
    xticks([])
    
    %% B3B6B9 Activity
    subplot(6,cols,i+6)
    
    % Model Data
    if i==2
        plot(T,out_lswallow.B6B9B3(i_ss_swallow),'Color',[72,129,191]/255,'LineWidth',LW)
        B6B9B3_model = interp1(T,out_lswallow.B6B9B3(i_ss_swallow),t_animal,"linear",0);
    else
        plot(T,out_uswallow.B6B9B3(i_ss_swallow),'Color',[72,129,191]/255,'LineWidth',LW)
        B6B9B3_model = interp1(T,out_uswallow.B6B9B3(i_ss_swallow),t_animal,"linear",0);
    end
    if i==1
        ylabel("B6/B9/B3","Rotation",rot,"HorizontalAlignment","right","VerticalAlignment","middle",'Color','k')
    end
    ylim([-0.1,1.1])
    yticks([0,1])
    
    % Animal Data
    yyaxis right
    hold all
    B6B9B3_animal = animal_data_loaded_swallow{4,1}{1,1};
    B6B9B3_animal_upper = animal_data_loaded_swallow{4,1}{2,1};
    B6B9B3_animal_lower = animal_data_loaded_swallow{4,1}{3,1};

    fill([t_animal,fliplr(t_animal)],[B6B9B3_animal_lower,fliplr(B6B9B3_animal_upper)],0.8*[72,129,191]/255,"EdgeColor","none","FaceAlpha",0.4)
    plot(t_animal,B6B9B3_animal,'-','Color',0.6*[72,129,191]/255,'LineWidth',LW)

    [r,lags] = xcorr(B6B9B3_animal,B6B9B3_model,"normalized");
    ind = find(lags==0);
    R_B6B9B3 = r(ind);

    yticks([])

    text(100,25,"R = "+num2str(R_B6B9B3,2))

    ylim([-0.1,1.1]*50)
    yticks([])
    
    xlim([0,100])
    xticks([])
    
    %% B4B5 Activity
    subplot(6,cols,i+8)
    
    % Model Data
    if i==2
        plot(T,out_lswallow.B4B5(i_ss_swallow),'Color',[8,171,111]/255,'LineWidth',LW)  
        B4B5_model = interp1(T,out_lswallow.B4B5(i_ss_swallow),t_animal,"linear",0);
    else
        plot(T,out_uswallow.B4B5(i_ss_swallow),'Color',[8,171,111]/255,'LineWidth',LW)  
        B4B5_model = interp1(T,out_uswallow.B4B5(i_ss_swallow),t_animal,"linear",0);
    end
    if i==1
        ylabel("B4/B5","Rotation",rot,"HorizontalAlignment","right","VerticalAlignment","middle",'Color','k')
    end
    ylim([-0.1,1.1])
    yticks([0,1])
    
    % Animal Data
    yyaxis right
    hold all
    B4B5_animal = animal_data_loaded_swallow{6,1}{1,1};
    B4B5_animal_upper = animal_data_loaded_swallow{6,1}{2,1};
    B4B5_animal_lower = animal_data_loaded_swallow{6,1}{3,1};

    fill([t_animal,fliplr(t_animal)],[B4B5_animal_lower,fliplr(B4B5_animal_upper)],0.8*[8,171,111]/255,"EdgeColor","none","FaceAlpha",0.4)
    plot(t_animal,B4B5_animal,'-','Color',0.6*[8,171,111]/255,'LineWidth',LW)

    [r,lags] = xcorr(B4B5_animal,B4B5_model,"normalized");
    ind = find(lags==0);
    R_B4B5 = r(ind);

    yticks([])

    text(100,12.5,"R = "+num2str(R_B4B5,2))

    ylim([-0.1,1.1]*20)
    yticks([])
    
    xlim([0,100])
    xticks([])
    
    %% Force Data
    subplot(6,cols,i+10)
    if i==2
        % Model Data
        F_model = (out_lswallow.F_fg(i_ss_swallow) + out_lswallow.F_fh(i_ss_swallow))/max((out_lswallow.F_fg(i_ss_swallow) + out_lswallow.F_fh(i_ss_swallow)));
        plot(T,F_model,'color',[0,0,0],'LineWidth',LW)
        
        
        % Animal Data
        % yyaxis right
        hold all
        F_animal = animal_data_loaded_swallow{7,1}{1,1};
        F_animal_upper = animal_data_loaded_swallow{7,1}{2,1};
        F_animal_lower = animal_data_loaded_swallow{7,1}{3,1};
        
        F_model = interp1(T,F_model,t_animal,"linear",0);
        [r,lags] = xcorr(F_animal,F_model,"normalized");
        ind = find(lags==0);
        R_F = r(ind);

        maxF = max(F_animal);
        minF = min(F_animal);
        fill([t_animal,fliplr(t_animal)],[(F_animal_lower - minF),fliplr((F_animal_upper - minF))]/(maxF - minF),0.6*[1,1,1],"EdgeColor","none","FaceAlpha",0.4)
        yticks([])
    
        text(100,.5,"R = "+num2str(R_F,2))
    end
    if i==1
        plot(T, 0*T,'-k','LineWidth',LW)
        ylabel("Force","Rotation",0,"HorizontalAlignment","right","VerticalAlignment","middle")
    end
    xlabel("Normalized Time [%Cycle]")
    ylim([-0.1,1.4])
    yticks([0,1])
    yticks([])
    xlim([0,100])
    
    for j=1:6
        subplot(6,cols,1+2*(j-1));
        set(gca,"FontSize",12)
        subplot(6,cols,2+2*(j-1));
        set(gca,"FontSize",12)
    end
end