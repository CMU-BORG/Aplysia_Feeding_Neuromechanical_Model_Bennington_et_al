%% Gill 2020 Loaded Swallow Data

file = "Gill2020_LoadedSwallows.csv";
folder = cd + "\Datasets\";

data = readtable(folder+file,"NumHeaderLines",1);

t = linspace(0,1,1000);

interp_data = cell(7,1);

lines = {'-k','--k','--k'};

t_0 = data{1,end-1};

data_temp = cell(3,1);
i=7; % force is the 7th dataset in this document
% figure; hold all
for j = 1:3
    t_ = data{:,6*(i-1) + 2*(j-1) + 1};
    x_ = data{:,6*(i-1) + 2*(j-1) + 2};
    
    nan_val = isnan(t_);
    t_(nan_val) = [];
    x_(nan_val) = [];

    [~,uniqueIdx] = unique(t_);
    t_ = t_(uniqueIdx);
    x_ = x_(uniqueIdx);
    x_([1,end]) = mean(x_([1,end])); % enforcing continuity between cycles

    t_ = [t_; t_+t_(end); t_+2*t_(end)];
    x_ = [x_;x_;x_];

    [~,uniqueIdx] = unique(t_);
    t_ = t_(uniqueIdx);
    x_ = x_(uniqueIdx);
        
    data_temp{j} = interp1(t_,x_,t+t_0,'');
%     plot(t+t_0,data_temp{j},lines{j})
%     plot(t_,x_,lines{j},'color','r')
    
end

med = data_temp{1};
q3 = data_temp{2};
q1 = data_temp{3};

% figure; hold all
% plot(t_,med(t_),'-k')
% plot(t_,q3(t_),'--k')
% plot(t_,q1(t_),'--k')



mean_val = (1/3)*(med + q1 + q3);
f_max =  max(mean_val);
mean_val = mean_val / f_max;

mean_force = @(t_) interp1(t,mean_val,t_,'spline');
std_force = @(t_) interp1(t,(1/1.35)*(q3 - q1)/f_max,t_,'spline');


