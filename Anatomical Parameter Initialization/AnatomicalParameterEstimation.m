%{
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Aplysia Californica Feeding Neuromechanical Model
                    Anatomical Parameter Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bennington, et al.

This script estimates anatomical parameters of the model from a high
resolution MRI image of the midsaggital plane reported in Neustadter et al
2007.
%}

% accessing functions from main folder
dir = split(cd,"\");
dir = strjoin(dir(1:end-1),"\");
addpath(dir)

%% Loading Image and Digitized Data
mri_img = im2double(imread('HighResolutionMidsagImage.PNG'));
mri_img = mri_img(:,:,1)*(1/3) + mri_img(:,:,2)*(1/3) + mri_img(:,:,3)*(1/3);

mri_data = readcell("HighResolutionMidsagData.csv");

% scale bar
scale_bar_x = [mri_data{3:4,1}];
scale_bar_y = [mri_data{3:4,2}];
px_to_cm = sqrt(diff(scale_bar_x).^2 + diff(scale_bar_y).^2);

% radular stalk width
RSW_x = [mri_data{3:4,7}];
RSW_y = [mri_data{3:4,8}];
RSW_px = sqrt(diff(RSW_x).^2 + diff(RSW_y).^2);

scale = RSW_px; % switch to px_to_cm to plot everything in cm
RSW = RSW_px / scale; 
V0 = 7.5*(RSW^3); % volume of the odontophore based on Neustadter et al 2002

[u,v] = size(mri_img);
R = imref2d(size(mri_img),[0,v/scale],[0,u/scale]); % used to scale image to same units as model


%% Rectangular I1I3 Approximation

% I1/I3 complex data
I1I3_x = [mri_data{3:6,3}] / scale;
I1I3_y = u/scale - [mri_data{3:6,4}] / scale;

x1_tilde = 0.5*(I1I3_x(1) + I1I3_x(2)); y1_tilde = 0.5*(I1I3_y(1) + I1I3_y(2));
x2_tilde = 0.5*(I1I3_x(2) + I1I3_x(3)); y2_tilde = 0.5*(I1I3_y(2) + I1I3_y(3));
x3_tilde = 0.5*(I1I3_x(3) + I1I3_x(4)); y3_tilde = 0.5*(I1I3_y(3) + I1I3_y(4));
x4_tilde = 0.5*(I1I3_x(4) + I1I3_x(1)); y4_tilde = 0.5*(I1I3_y(4) + I1I3_y(1));

% Rest dimensions of the I1/I3 lumen
H_lumen = y2_tilde - y4_tilde;
L_lumen = x1_tilde - 0.5*(I1I3_x(3) + I1I3_x(4));

H_lumen_tracing = (0.6 - 0.3177); % [units] 
tracing_scale_factor = H_lumen / H_lumen_tracing; % [R / unit]

L_head = 0.82*tracing_scale_factor; % [R] based on head outline from Webster-Wood et al. 2020


% initial position of the lateral groove end point w.r.t. jaw line
xd_tilde = x1_tilde - I1I3_x(3);
xv_tilde = x1_tilde - I1I3_x(4); 

lumen = polyshape([x1_tilde,x1_tilde,I1I3_x(3),I1I3_x(4)], ...
                  [y4_tilde,y4_tilde+H_lumen,y4_tilde+H_lumen,y4_tilde] - y4_tilde);


%% Odontophore Fitting

% odontophore outline
odont_x = [mri_data{3:end,5}] / scale;
odont_y = u/scale - [mri_data{3:end,6}] / scale - y4_tilde;

idx_nan = isnan(odont_x);
odont_x(idx_nan) = []; odont_y(idx_nan) = [];

% creating a callable function that gives the value of the odontophore data
% at a given angle
xg_data = mean(odont_x); yg_data = mean(odont_y);
theta_data = atan2(odont_y - yg_data,odont_x - xg_data);
R_data = sqrt((odont_x - xg_data).^2 + (odont_y - yg_data).^2);
theta_data = [theta_data(1:end-1)-2*pi, theta_data, theta_data(2:end)+2*pi];
[theta_data,u_idx] = unique(theta_data);
R_data = [R_data(1:end-1), R_data, R_data(2:end)];
R_data = R_data(u_idx);
R_data_func = @(x) interp1(theta_data,R_data,x,'spline');
X_data_func = @(x) [xg_data;yg_data] + R_data_func(x).*[cos(x);sin(x)];

% minimal kinematic parameters to calculate volume
kinematicParams = struct();
left_curve = linspace(0,0.05,5);
right_curve = linspace(0.95,1,5);
mid_curve = linspace(0.05,0.95,22);
mid_curve([1,end]) = [];

kinematicParams.M = 30;
kinematicParams.N = 15;
kinematicParams.s_vec = [left_curve,mid_curve,right_curve];
kinematicParams.V0 = V0;

% initial guess for odontophore parameters based on data
xg0 = xg_data; yg0 = yg_data;
x_std = std(odont_x - xg_data); y_std = std(odont_y - yg_data);
L0 = 0.75*sqrt(x_std.^2 + y_std.^2);
thetag0 = atan2(y_std,x_std)*1.25;
R10 = mean([x_std,y_std]); R20 = R10;

kinematicParams.phi = GetTangent(R10,R20,L0);

% function enforcing the volume constraint to obtain alpha
f_vol = @(alpha) V0 - OdontophoreVolume([R10,R20,L0,alpha],kinematicParams);
alpha = fsolve(f_vol,1,optimset('Display','off'));

% error function to fit odontophore model to high res MRI 
f_err = @(X) OdontophoreError_com(X_data_func,X(1),X(2),X(3),X(4),X(5),X(6),V0,kinematicParams);
X = fminsearch(f_err,[R10,R20,L0,xg0,yg0,thetag0],optimset("TolFun",1e-6,"TolX",1e-6));

% Initial odontophore parameters for dynamical model
R10 = X(1); R20 = X(2); L0 = X(3);
xg0 = X(4); yg0 = X(5); thetag = X(6);
xg0_tilde = x1_tilde - X(4);
yg0_tilde = yg0; thetag0_tilde = thetag;
f_vol = @(alpha) V0 - OdontophoreVolume([R10,R20,L0,alpha],kinematicParams);
alpha0 = fsolve(f_vol,1,optimset('Display','off'));

[x_pts,y_pts] = DrawOdontophoreMS([xg0,yg0,thetag0,R10,R20,L0,alpha0],kinematicParams);
odont = polyshape(x_pts,y_pts);

% finding the relative position of the hinge on the odontophore
s_hinge_func = @(s) Pos_Global(s,1,-pi/2,[xg0,yg0,thetag0,R10,R20,L0,alpha0],kinematicParams,1)'*[0;1;0];
s_hinge = fsolve(s_hinge_func,0.24,optimset("TolFun",1e-15,"MaxFunEvals",1e8,"MaxIter",1e8));

xyz_hinge = Pos_Global(s_hinge,1,-pi/2,[xg0,yg0,thetag0,R10,R20,L0,alpha0],kinematicParams,1);


%% Head Data
% reading in data points for the head outline
file = 'HeadOutlineData.csv';
data = readmatrix(file);
data(:,1) = tracing_scale_factor*data(:,1);
data(:,2) = tracing_scale_factor*0.75*data(:,2);

jaw_v_line = 2.5; % position of the ventral jaw line in the scaled tracing
data(:,1) = data(:,1) - L_head + x1_tilde;
data(:,2) = data(:,2) - jaw_v_line; % registering outline so ventral jaw line is at y=0

head_x = [data(:,1); data(1,1)];
head_y = [data(:,2); data(1,2)];

%% Plotting the Fitted Model on top of the High Resolution MRI
figure('Color','w')
% MRI Image
imshow(flipud(mri_img),'XData',[0,v/scale],'YData',[0,u/scale] - y4_tilde); 
set(gca, 'YDir', 'normal');
set(gca,'FontName','Arial')
axis on;
hold all
xlabel("X Position [RSW]")
ylabel("Y Position [RSW]")

% Head outline
% plot(head_x,head_y,'-m','LineWidth',2,'HandleVisibility','off')

% Model Approximation of Rest Geometry
plot(lumen,'EdgeColor','k','FaceColor',[75,0,130]/255,'DisplayName',"Model I1/I3 Lumen")
plot(odont,'EdgeColor','k','FaceColor',[1,0,0.6],'DisplayName',"Model Odontophore")
plot(xyz_hinge(1),xyz_hinge(2),'sk','MarkerFaceColor','r','DisplayName','Model Hinge Point')
plot(xg0,yg0,'sk','MarkerFaceColor','k','DisplayName','Model Center of Mass')

% Digitize Data
% plot(scale_bar_x,scale_bar_y,'or','MarkerFaceColor','r','HandleVisibility','off')
plot(I1I3_x,I1I3_y - y4_tilde,'ok','MarkerFaceColor',[130,0,200]/255,'DisplayName',"I1/I3 Complex Bound"+newline+"Box Corners");
plot(odont_x,odont_y,'ok','MarkerFaceColor','m','DisplayName',"Odontophore Outline")
plot(RSW_x / scale,u/scale - RSW_y / scale  - y4_tilde,'ok','MarkerFaceColor','b','DisplayName','RSW Reference Points')

legend('NumColumns',4)

%% Saving Anatomical Parameters to a struct for later use
anatomical_parameters = struct();

% reference geometry values
anatomical_parameters.L_head = L_head;
anatomical_parameters.L_lumen = L_lumen;
anatomical_parameters.H_lumen = H_lumen;
anatomical_parameters.s_hinge = s_hinge;
anatomical_parameters.V0 = V0;

% parameters for state initialization
anatomical_parameters.xv_tilde = xv_tilde;
anatomical_parameters.xd_tilde = xd_tilde;
anatomical_parameters.xg0_tilde = xg0_tilde;
anatomical_parameters.yg0_tilde = yg0_tilde;
anatomical_parameters.thetag0_tilde = thetag0;
anatomical_parameters.R10 = R10;
anatomical_parameters.R20 = R20;
anatomical_parameters.L0 = L0;
anatomical_parameters.alpha0 = alpha0;

anatomical_parameters.head_x = head_x;% - x1_tilde;
anatomical_parameters.head_y = head_y;

save("AnatomicalParams","anatomical_parameters");
