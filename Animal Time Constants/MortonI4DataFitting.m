%% I4 Time Constant Fitting
% clear all; close all; clc;
% dir = split(cd,"\");
% dir = strjoin(dir(1:end-1),"\");
% addpath(dir+"\Animal Data\Datasets\");
% 
data = readmatrix("MortonI4Data.csv","NumHeaderLines",2);
%%
t1 = data(:,1);
F1 = data(:,2);
idx = isnan(t1);
t1(idx) = []; F1(idx) = [];
t1 = t1 - t1(1); F1 = F1 - F1(1); F1 = F1 / max(F1);

t2 = data(:,3);
F2 = data(:,4);
idx = isnan(t2);
t2(idx) = []; F2(idx) = [];
t2 = t2 - t2(1); F2 = F2 - F2(1); F2 = F2 / max(F2);

t3 = data(:,5);
F3 = data(:,6);
idx = isnan(t3);
t3(idx) = []; F3(idx) = [];
t3 = t3 - t3(1); F3 = F3 - F3(1); F3 = F3 / max(F3);

pulses = data(1:6,7);
dT = 4; %diff(pulses); dT = mean(dT([1,3,5])); 

figure("Position",[100,100,450,400],"Color","w"); hold all
scatter(t1,F1,'o','MarkerEdgeColor','none','MarkerFaceColor','r','DisplayName','Cycle 1','MarkerFaceAlpha',0.4);
scatter(t2,F2,'o','MarkerEdgeColor','none','MarkerFaceColor','g','DisplayName','Cycle 2','MarkerFaceAlpha',0.4);
scatter(t3,F3,'o','MarkerEdgeColor','none','MarkerFaceColor','b','DisplayName','Cycle 3','MarkerFaceAlpha',0.4);

t_all = [t1;t2;t3]; 
F_all = [F1;F2;F3];
tt = linspace(0,max(t_all),200);
ta = 0.2;
% plot(tt,(tt>=ta).*(tt<=ta+dT),'--','Color',0.6*[1,1,1])


err = @(X) sum( (F_all - Tension(t_all,X(1),X(2),X(2)+dT)).^2 );
X = fminsearch(err,[0.5,0.5,1]);

T_ = tt;
% ta = 1.75; dt = 0.2; tb = ta + dt;
% tau1 = 0.2; tau2 = 0.6;
ta = X(2); tb = X(2) + dT;
tau = X(1);

plot(T_,Tension(T_,tau,ta,tb),'-','LineWidth',2,"DisplayName","Model",'Color',0.*[1,1,1]);
fill([ta,ta,tb,tb],[0,1,1,0],0.6*[1,1,1],"FaceAlpha",0.2,"EdgeColor","none","DisplayName","Stimulation On")

legend("Location","best");
xlabel("Time [s]")
ylabel("Normalized Response [ ]")
yticks([0:0.2:1])
title("I4 - Morton et al. 1993")
ylim([0,1])
xlim([0,max(T_)])
set(gca,"FontName","Arial","FontSize",12)


fprintf("I4 data:\n")
fprintf("    tau [s] = %.3f\n",tau)





function T = Tension(t,tau,t_on,t_off)


T_rise = @(t) ( 1 - exp(-(t-t_on)/tau) .* (1 + (t-t_on)/tau) );
T_peak = T_rise(t_off);
A_peak = 1 - exp(-(t_off-t_on)/tau);

T_fall = @(t) exp(-(t-t_off)/tau) .* (T_peak + A_peak*(t-t_off)/tau);

T = (t>t_on).*(t<=t_off).*T_rise(t) + (t>t_off).*T_fall(t);
T = T / max(T);

end
