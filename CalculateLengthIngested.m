%% Calculate Length of Ingested Seaweed
%{
    Calculation of length of ingested seaweed that smooths jitteriness in
    still solution.
%}

ds = out; % dataset to use for calculation

t_ = ds.tout;   % time vector of solution

xg = ds.xg;     % center of mass of odontophore
yg = ds.yg;
xh = ds.xh;     % head position
thetag = ds.theta_g;                        % angle of the odontophore

xf = -(L_head + xh) + (xg - grasperR*cos(thetag - theta_fg));  % position of the food attachment point relative to the head

fg_state = ds.fg_state;                                         % grasper friction state 
FI3_ant_nom = ds.P_I3ant;                                       % normalized force in the I3 anterior region (to tell if jaws are closed)
out_sens_mechanical_grasper = ds.out_sens_mechanical_grasper;   % presence of food at the grasper
fixation = ds.fixation;

sxf = interp1(t_,smooth(xf,100),"spline","pp"); % spline interpolate the food position to remove infinitessimal jitters
dxf = fnder(sxf,1);                             % first derivative of the spline interpolation
dxf = ppval(dxf,t_);                            % first derivative calculated at all time points

grasper_closed_on_food = fg_state .* out_sens_mechanical_grasper .* (1-fixation); % is the grasper closed onto food?
pushing_forward = dxf > 0; % is the food moving forward?
jaws_closed = FI3_ant_nom > 0.5; % are the jaws closed to prevent egestion?

% indicator function for whether or not to integrate motion (called
% delta_ingest(tau) in the Supplemental Information)
food_moving = pushing_forward.*( grasper_closed_on_food .* (1-jaws_closed) ) + (1 - pushing_forward).*( grasper_closed_on_food );

% calculating the time integral of seaweed ingestion
dL_ingested = food_moving .* dxf;
L_ingested = cumtrapz(t_,-dL_ingested);

% saving to struct
out.L_ingested = L_ingested;

