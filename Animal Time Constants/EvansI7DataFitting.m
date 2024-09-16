%% I7 Time Constant Fitting
% clear all; close all; clc;
% dir = split(cd,"\");
% dir = strjoin(dir(1:end-1),"\");
% addpath(dir+"\Animal Data\Datasets\");

data = readmatrix("EvansI7Data.csv","NumHeaderLines",2);

%%

t = data(:,1);
T = data(:,2);
spikes = [data(1:5,3);max(t)];

min_t = 0;
min_t_end = 1000;
% figure; hold all

Ts = cell(5,1);
for i=1:5
    idx = boolean((t>=spikes(i)).*(t<=spikes(i+1)));
    t_ = t(idx);
    T_ = T(idx);
    
    pk = find(T_==max(T_),1);
    t_ = t_ - t_(pk);

    if min(t_)<min_t
        min_t = min(t_);
    end

    if max(t_)<min_t_end
        min_t_end = max(t_);
    end

    T_ = T_ - T_(1);
    T_ = T_ / max(T_);

    Ts{i} = [t_, T_];

%     plot(t_,T_,'.')
end
%%
colors = {'c','m','g','b','r'};
figure; hold all
tt = linspace(0,min_t_end,200);
t_all = [];
T_all = [];
for i=1:5
    data = Ts{i};
    t_ = data(:,1) - min_t;
    idx = boolean(t_>min_t_end - min_t);
    t_(idx) = [];
    t_ = [linspace(min_t,t_(1),50)';t_];
    [t_,ui] = unique(t_);
    
    T_ = data(:,2);
    T_(idx) = [];
    T_ = [zeros(50,1);T_];
    T_ = T_(ui);

    TT = interp1(t_,T_,tt,'pchip');
    scatter(tt,TT,'o','MarkerEdgeColor','none','MarkerFaceColor',colors{i},'DisplayName',"Cycle "+num2str(i),"MarkerFaceAlpha",0.4)

    t_all = [t_all,tt];
    T_all = [T_all,TT];
    
end

err = @(X) sum( (T_all - Tension(t_all,X(1),X(2),X(3),X(3)+0.001)).^2 );
X = fminsearch(err,[0.2,0.6,1.75,1]);

T_ = tt;
% ta = 1.75; dt = 0.2; tb = ta + dt;
% tau1 = 0.2; tau2 = 0.6;
ta = X(3); tb = X(3) + 0.001;
tau1 = X(1); tau2 = X(2);

plot(T_,Tension(T_,tau1,tau2,ta,tb),'-','LineWidth',2,"DisplayName","Model",'Color',0.*[1,1,1]);


legend("Location","best");
xlabel("Time [s]")
ylabel("Normalized Response [ ]")
yticks([0:0.2:1])
% title("Collostylar Cap - Evans et al. 1996")

fprintf("Collostylar Cap data:\n")
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