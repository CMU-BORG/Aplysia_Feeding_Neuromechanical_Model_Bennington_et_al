%% Fitting time constants for hinge from Sutton et al

% clear all; close all; clc;
% dir = split(cd,"\");
% dir = strjoin(dir(1:end-1),"\");
% addpath(dir+"\Animal Data\Datasets\");
data = readmatrix("SuttonHingeData.csv","NumHeaderLines",2);

t0_1r = data(1,1);
t0_2r = data(3,1);
t0_1f = data(2,1);
t0_2f = data(4,1);
dT = 0.5*( (t0_1f - t0_1r) + (t0_1f - t0_1r) ); % duration of stimulation


t_1r = data(:,3);
idx1 = isnan(t_1r);
t_1r(idx1) = [];
F_1r = data(:,4);
F_1r(idx1) = [];

t_1r = t_1r - t0_1r;
F_1r = F_1r - F_1r(1);
F_1r = F_1r / F_1r(end);

t_2r = data(:,5);
idx2 = isnan(t_2r);
t_2r(idx2) = [];
F_2r = data(:,6);
F_2r(idx2) = [];

t_2r = t_2r - t0_2r;
F_2r = F_2r - F_2r(1);
F_2r = F_2r / F_2r(end);





t_1f = data(:,7);
idx1 = isnan(t_1f);
t_1f(idx1) = [];
F_1f = data(:,8);
F_1f(idx1) = [];

t_1f = t_1f - t0_1r;
F_1f = F_1f - F_1f(end);
F_1f = F_1f / F_1f(1);

t_2f = data(:,9);
idx2 = isnan(t_2f);
t_2f(idx2) = [];
F_2f = data(:,10);
F_2f(idx2) = [];

t_2f = t_2f - t0_2r;
F_2f = F_2f - F_2f(end);
F_2f = F_2f / F_2f(1);


T1 = [t_1r;t_1f];
T2 = [t_2r;t_2f];
F1 = [F_1r;F_1f];
F2 = [F_2r;F_2f];

idx1 = (T1<20);
idx2 = (T2<20); % getting rid of excessive time points
T1 = T1(idx1);
F1 = F1(idx1);
T2 = T2(idx2);
F2 = F2(idx2);
T = [T1;T2];
F = [F1;F2];

figure("Position",[100,100,450,400],"Color","w"); hold all
scatter(T1,F1,'o','MarkerEdgeColor','none','MarkerFaceColor','g','DisplayName',"Cycle "+num2str(1),"MarkerFaceAlpha",0.4)
scatter(T2,F2,'o','MarkerEdgeColor','none','MarkerFaceColor','b','DisplayName',"Cycle "+num2str(2),"MarkerFaceAlpha",0.4)


err = @(X) sum( (F - Tension(T,X(1),X(2),X(2)+dT)).^2 );
X = fminsearch(err,[.1,0]);
% X(2) = 0;

T_ = linspace(0,max(T),200);
ta = X(2); tb = X(2) + dT;
tau = X(1); 
plot(T_,Tension(T_,tau,ta,tb),'-','LineWidth',2,"DisplayName","Model",'Color',0.*[1,1,1]);
fill([ta,ta,tb,tb],[0,1,1,0],0.6*[1,1,1],"FaceAlpha",0.2,"EdgeColor","none","DisplayName","Stimulation On")

legend("Location","best");
xlabel("Time [s]")
ylabel("Normalized Response [ ]")
title("Hinge - Sutton et al. 2004")
yticks([0:0.2:1])
ylim([0,1])
xlim([0,max(T)])
set(gca,"FontName","Arial","FontSize",12)

fprintf("Hinge data:\n")
fprintf("    tau [s] = %.3f\n",tau)

function T = Tension(t,tau,t_on,t_off)


T_rise = @(t) ( 1 - exp(-(t-t_on)/tau) .* (1 + (t-t_on)/tau) );
T_peak = T_rise(t_off);
A_peak = 1 - exp(-(t_off-t_on)/tau);

T_fall = @(t) exp(-(t-t_off)/tau) .* (T_peak + A_peak*(t-t_off)/tau);

T = (t>t_on).*(t<=t_off).*T_rise(t) + (t>t_off).*T_fall(t);
T = T / max(T);

end