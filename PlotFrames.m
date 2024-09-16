%% Plotting Frames of the model 
%{
    Generate and plot the three configurations of the model shown in Fig 2c
%}

data = out_bite;

N_frames = 3;

%% Get 3 equally space frames between peak retraction and peak protraction
x_gh = data.x_gh(i_ss_bite);
[~,min_xgh] = min(x_gh);
[~,max_xgh] = max(x_gh);

t_range = data.tout(i_ss_bite(min_xgh:max_xgh));
t_frames = linspace(t_range(1),t_range(end),N_frames);

%% Plotting Model
f1 = figure("Position",[100,100,1500,300],"Color","w");

for i=1:N_frames
    subplot(1,N_frames,i); hold all
    [~,t_ind] = min( abs(data.tout(i_ss_bite) - t_frames(i) ) );
    ind = i_ss_bite(t_ind);

    Xi = [data.xg(ind);
          data.yg(ind);
          data.theta_g(ind);
          data.xh(ind);
          data.xd(ind);
          data.xv(ind)];
    Ai = 0.5*ones(9,1); Ai(4) = 0.2;

    I2_attach = [data.Xtd(ind,:),data.Xtv(ind,:)];

    f1 = PlotModel(gca,Xi,head_outline,buccal_mass_params,I2_attach,Ai);
    ylim([-2.5,6.5])

end