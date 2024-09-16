%% Fitting time constants for anterior I3 pinch from McManus et al 2014

% clear all; close all; clc;


data = readmatrix("McManusI3aData.csv");

%%
t = data(:,1);
f = data(:,2);

T = [];
F = [];
[~,cycle_starts] = findpeaks(-f);
cycle_starts = [1;cycle_starts;length(t)];

colors = ['r','g','b'];

figure("Position",[100,100,450,400],"Color","w"); hold all
for i=1:3
    idx = cycle_starts(i)+1*(i==3):cycle_starts(i+1); % removing an extra point from 3rd cycle
    
    t_ = t(idx); f_ = f(idx);
    t_ = t_ - t_(1);
    f_ = f_-f_(1); f_ = f_ / max(f_);

    T = [T;t_];
    F = [F;f_];
    
    scatter(t_,f_,'o','MarkerEdgeColor','none','MarkerFaceColor',colors(i),'DisplayName',"Cycle "+num2str(i),"MarkerFaceAlpha",0.4)
end

t_on = 0;

% 2 seconds because that was used in the protocol
err = @(X) sum( (F - Tension(T,X(1),X(2),X(2)+2)).^2 ) + (1e10)*(X(2) < 0);
X = fminsearch(err,[.1,2]);
% X(2) = 0;5

T_ = linspace(0,max(T),200);
ta = X(2); tb = X(2) + 2;
tau = X(1); 
plot(T_,Tension(T_,tau,ta,tb),'-','LineWidth',2,"DisplayName","Model",'Color',0.*[1,1,1]);
fill([ta,ta,tb,tb],[0,1,1,0],0.6*[1,1,1],"FaceAlpha",0.2,"EdgeColor","none","DisplayName","Stimulation On")

legend("Location","best");
xlabel("Time [s]")
ylabel("Normalized Response [ ]")
title("I3 pinch - McManus et al. 2014")
yticks([0:0.2:1])
ylim([-0.,1.])
xlim([0,max(T)])
set(gca,"FontName","Arial","FontSize",12)


fprintf("I3 pinch data:\n")
fprintf("    tau [s] = %.3f\n",tau)

function T = Tension(t,tau,t_on,t_off)


T_rise = @(t) ( 1 - exp(-(t-t_on)/tau) .* (1 + (t-t_on)/tau) );
T_peak = T_rise(t_off);
A_peak = 1 - exp(-(t_off-t_on)/tau);

T_fall = @(t) exp(-(t-t_off)/tau) .* (T_peak + A_peak*(t-t_off)/tau);

T = (t>t_on).*(t<=t_off).*T_rise(t) + (t>t_off).*T_fall(t);
T = T / max(T);

end