%% SRF Time Constant Fitting
% clear all; close all; clc;
% dir = split(cd,"\");
% dir = strjoin(dir(1:end-1),"\");
% addpath(dir+"\Animal Data\Datasets\");

data = readmatrix("BorovikovSRFData.csv");

%%

t = data(:,1); t = t - t(1);
T = data(:,2);
dT = 1; % how long the stimulus was


figure; hold all
scatter(t,T,'o',"MarkerEdgeColor","none",'MarkerFaceColor','r',"DisplayName","Cycle 1",'MarkerFaceAlpha',0.4)


err = @(X) sum( (T - Tension(t,X(1),X(2),X(3),X(3)+dT)).^2 );
X = fminsearch(err,[0.5,1,0]);

T_ = linspace(0,max(t),200);
ta = X(3); tb = X(3) + dT;
tau1 = X(1); tau2 = X(2);

plot(T_,Tension(T_,tau1,tau2,ta,tb),'-','LineWidth',2,"DisplayName","Model",'Color',0.*[1,1,1]);

legend("Location","best");
xlabel("Time [s]")
ylabel("Normalized Response [ ]")
yticks([0:0.2:1])
% title("Subradular Fiber Data - Borovikov 2000")

fprintf("SRF data:\n")
fprintf("    tau [s] = %.3f\n",tau1)
fprintf("    beta [ ] = %.3f\n",tau1 / tau2)



function T = Tension(t,tau1,tau2,t_on,t_off)


T_rise = @(t) ( 1 - exp(-(t-t_on)/tau1) .* (1 + (t-t_on)/tau1) );
T_peak = T_rise(t_off);
A_peak = 1 - exp(-(t_off-t_on)/tau1);

T_fall = @(t) exp(-(t-t_off)/tau2) .* (T_peak + A_peak*(t-t_off)/tau2);

T = (t>t_on).*(t<=t_off).*T_rise(t) + (t>t_off).*T_fall(t);
T = T / max(T);

end