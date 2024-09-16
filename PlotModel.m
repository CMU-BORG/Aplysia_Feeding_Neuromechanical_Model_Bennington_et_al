function fig = PlotModel(f1,X_state, outline_data, buccal_mass_data, I2_attach, muscle_activations)
%% Function to plot the feeding Aplysia model given the current state
%{
    X_state:            state variables of the model (x_g,y_g,theta_g,x_h,x_d,x_v)
    outline_data:       matrix to plot the Aplysia head outline
    buccal_mass_data:   design parameters of the buccal mass:
                        - odontophore radius (R)
                        - relative (to head) jaw x pos (jaws_x)
                        - ventral jaw y pos (jaws_v_y)
                        - dorsal jaw y pos (jaws_d_y)
                        - E6 anchor x
                        - theta_H angle from rad stalk to hinge
                        - theta_fg angle from rad stalk to food attachment
    I2_attach:          vector of I2 odontophore attachment points
                        - [xd_I2, yd_I2, xv_I2, yv_I2]
    muscle_activations: vector of normalized muscle activation levels
                        - I1d, I1v, I2, I3, Hinge, E1, E6
%}

%% Meta parameters to tweak
ylims = [-0.2,4];

%% Separating input variable vectors 

% Separating out the state variables
x_g = X_state(1);
y_g = X_state(2);
theta_g = X_state(3);
x_h = X_state(4);
x_d = X_state(5);
x_v = X_state(6);

% separating the head outline data matrix
head_x = [outline_data(:,1); outline_data(1,1)]; 
head_y = [outline_data(:,2); outline_data(1,2)];

% separating out the buccal mass data 
R = buccal_mass_data(1);            % radius of odontophore
jaws_x = buccal_mass_data(2);       % relative x pos of jaw line
jaws_v_y = buccal_mass_data(3);     % relative y pos of ventral jaw line
jaws_d_y = buccal_mass_data(4);     % relative y pos of dorsal jaw line
E6_anchor_x = buccal_mass_data(5);  % relative E6 anchor position
theta_H = buccal_mass_data(6);
theta_fg = buccal_mass_data(7);

del_x = linspace(-R,R,100);             % x points of the circle centered at origin (for odontophore
lumen_mid = mean([jaws_v_y,jaws_d_y]);  % height of the lumen mid line

I1d_act = muscle_activations(1);
I1v_act = muscle_activations(2);
I2_act = muscle_activations(3);
I3_act = muscle_activations(4);
I3ant_act = muscle_activations(5);
H_act = muscle_activations(6);
E1_act = muscle_activations(7);
E2_act = muscle_activations(8);
E6_act = muscle_activations(9);

%% Defining the visualization variables

def_thick = 1.5;
line_thick = @(a) 0.6*def_thick + def_thick*a; % defining thickness of line element based on muscle activation
line_alpha = @(a) 0.6 + 0.4*a;

I1_color = 0.7*[90/255, 155/255, 197/255];
I2_color = [220/255, 81/255, 81/255];
I3_color = [90/255, 155/255, 197/255];
H_color = [70/255, 185/255, 255/255];
E1_color = [125/255, 170/255, 50/255];
E2_color = [180/255, 230/255, 120/255];
E6_color = [50/255, 170/255, 125/255];

jaw_color = 0.8*[246/255, 253/255, 89/255];
ref_color = [100/255, 100/255, 100/255];

head_color = [0.9 0.6 0.6];

%% Setting up figure for plotting
% fig = figure(f1);
fig = -1;
cla reset
hold all; axis equal
ylim(ylims);

%% Plotting the head outline and springs
plot(x_h + head_x,head_y,'-','Color',head_color,'LineWidth',def_thick,'HandleVisibility','off')                                                            % head outline

spring_x = mean([head_x(1),head_x(end-1)]);
spring_y = mean([head_y(1),head_y(end-1)]);
plot([spring_x,x_h + spring_x],[spring_y,spring_y],'-o','Color',head_color,'MarkerFaceColor',head_color,'HandleVisibility','off')

%% Plotting the jaw points and reference lines (fixed)
plot(x_h + jaws_x,jaws_d_y,'o','Color',jaw_color,'MarkerFaceColor',jaw_color,'HandleVisibility','off');        % dorsal jaw anchor

plot(x_h + jaws_x,jaws_v_y,'o','Color',jaw_color,'MarkerFaceColor',jaw_color,'HandleVisibility','off');        % ventral jaw anchor

yline(jaws_d_y,'--','Color',ref_color,'LineWidth',0.5*def_thick,'HandleVisibility','off')                                           % top of lumen
yline(jaws_v_y,'--','Color',ref_color,'LineWidth',0.5*def_thick,'HandleVisibility','off')                                           % bottom of lumen and pivot rail
yline(lumen_mid,'--','Color',ref_color,'LineWidth',0.5*def_thick,'HandleVisibility','off')                                          % lumen midline (line of action for I3)

l1 = plot([x_h + jaws_x,x_h + jaws_x],[jaws_d_y,jaws_v_y],'--','Color',jaw_color,'LineWidth',line_thick(I3ant_act),'DisplayName','I3 Pinch'); % jaw line
l1.Color(4) = line_alpha(0.5*I3ant_act);

%% Plotting the I1 
plot(x_d,jaws_d_y,'o','Color',I1_color,'MarkerFaceColor',I1_color,'HandleVisibility','off')                               % dorsal I1 anchor
plot(x_v,jaws_v_y,'o','Color',I1_color,'MarkerFaceColor',I1_color,'HandleVisibility','off')                               % ventral I1 anchor

l1 = plot([x_d, x_h+jaws_x],[jaws_d_y,jaws_d_y],'-','Color',I1_color,'LineWidth',line_thick(I1d_act),'DisplayName','Dorsal I1');
l1.Color(4) = line_alpha(I1d_act);

l1 = plot([x_v, x_h+jaws_x],[jaws_v_y,jaws_v_y],'-','Color',I1_color,'LineWidth',line_thick(I1v_act),'DisplayName','Ventral I1');
l1.Color(4) = line_alpha(I1v_act);

% lateral groove
plot([x_d,x_v],[jaws_d_y,jaws_v_y],'--','Color',I1_color,'LineWidth',def_thick,'HandleVisibility','off')

%% Plotting the I3
l1 = patch([jaws_x + x_h,jaws_x + x_h,x_d,x_v],[jaws_v_y,jaws_d_y,jaws_d_y,jaws_v_y],I3_color,'DisplayName','I3');
l1.FaceAlpha = line_alpha(0.1*I3_act);
l1.EdgeColor = 'none';
%% Plotting the odontophore
x_circ = x_g-del_x;                             % x points of the odontophore
y_circ = real(sqrt(R^2 - (x_circ-x_g).^2));     % y points of the odontophore
plot(x_g,y_g,'o','Color',head_color,'LineWidth',def_thick,'HandleVisibility','off')        % center of mass of odontophore
plot(x_circ,y_g + y_circ,'-','Color',head_color,'LineWidth',def_thick,'HandleVisibility','off')    % plotting odontophore
plot(x_circ,y_g - y_circ,'-','Color',head_color,'LineWidth',def_thick,'HandleVisibility','off')

%% odontophore midline axis (to define the angle of the odontophore)
rad_stalk_axis_x = [x_g , x_g + R*cos(theta_g)];
rad_stalk_axis_y = [y_g , y_g + R*sin(theta_g)];
plot(rad_stalk_axis_x,rad_stalk_axis_y,'--','Color',ref_color,'LineWidth',0.5*def_thick,'HandleVisibility','off')

%% plotting food attachment
food_axis_x = [x_g , x_g + R*cos(theta_g - theta_fg)];
food_axis_y = [y_g , y_g + R*sin(theta_g - theta_fg)];
plot(food_axis_x,food_axis_y,'-o','Color',ref_color,'MarkerFaceColor',ref_color,'LineWidth',0.5*def_thick,'HandleVisibility','off')


%% plotting the I2

% dorsal tangent point
Xtd = [I2_attach(1), I2_attach(2)];
theta1 = acos((Xtd(1) - x_g)/R);
if ~isreal(theta1)
    disp('Theta 1 is complex.')
end

plot(Xtd(1),Xtd(2),'x','Color','k','MarkerFaceColor','k','LineWidth',1.5,'HandleVisibility','off')
l1 = line([Xtd(1),x_d],[Xtd(2),jaws_d_y],'Color',I2_color,'LineWidth',line_thick(I2_act),'HandleVisibility','off');
l1.Color(4) = line_alpha(I2_act);

% ventral tangent point
Xtv = [I2_attach(3), I2_attach(4)];
theta2 = (pi + acos((x_g - Xtv(1))/R));

plot(Xtv(1),Xtv(2),'x','Color','k','MarkerFaceColor','k','LineWidth',1.5,'HandleVisibility','off')
l1 = line([Xtv(1),x_v],[Xtv(2),jaws_v_y],'Color',I2_color,'LineWidth',line_thick(I2_act),'HandleVisibility','off');
l1.Color(4) = line_alpha(I2_act);

% wrapping around odontophore
theta_range = linspace(theta1,theta2,75);
I2_x = x_g + R*cos(theta_range);
I2_y = y_g + R*sin(theta_range);

l1 = plot(I2_x,I2_y,'-','Color',I2_color,'MarkerFaceColor',I2_color,'LineWidth',line_thick(I2_act),'DisplayName','I2');
l1.Color(4) = line_alpha(I2_act);



%% odontophore hinge attachment
x_hinge = x_g + R*cos(theta_g + theta_H); 
y_hinge = y_g + R*sin(theta_g + theta_H);
plot(x_hinge,y_hinge,'o','Color',H_color,'MarkerFaceColor',H_color,'HandleVisibility','off')
plot([x_g,x_hinge],[y_g,y_hinge],'--','Color',H_color,'HandleVisibility','off')

l1 = plot([x_hinge, x_v],[y_hinge, jaws_v_y],'-','Color',H_color,'LineWidth',line_thick(H_act),'DisplayName','Hinge');
l1.Color(4) = line_alpha(H_act);


%% E1 
l1 = plot([Xtd(1),x_h + jaws_x],[Xtd(2),jaws_d_y],'-','Color',E1_color,'LineWidth',line_thick(E1_act),'DisplayName','E1');
l1.Color(4) = line_alpha(E1_act);

%% E2 
l1 = plot([Xtv(1),x_h + jaws_x],[Xtv(2),jaws_v_y],'-','Color',E2_color,'LineWidth',line_thick(E2_act),'DisplayName','E2');
l1.Color(4) = line_alpha(E2_act);

%% E6
lat_groove_x = mean([x_d,x_v]);
lat_groove_y = mean([jaws_d_y,jaws_v_y]);

l1 = plot([lat_groove_x,x_h+E6_anchor_x],[lat_groove_y,jaws_d_y],'-','Color',E6_color,'LineWidth',line_thick(E6_act),'DisplayName','E6');
l1.Color(4) = line_alpha(E6_act);

end
