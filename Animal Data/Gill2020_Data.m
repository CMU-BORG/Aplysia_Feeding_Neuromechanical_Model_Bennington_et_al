%% Gill 2020 Loaded Swallow Data

loaded = 0;

if loaded
    file = "Gill2020_LoadedSwallows.csv";
    outfile = "Gill2020LoadedSwallow.mat";
    max_ind = 7;
else
    file = "Gill2020_UnloadedSwallows.csv";
    outfile = "Gill2020UnloadedSwallow.mat";
    max_ind = 6;
end
folder = cd + "\Datasets\";

data = readtable(folder+file,"NumHeaderLines",1);

t = linspace(0,1,1000);

lines = {'-k','--k','--k'};
t_min = Inf; t_max = -Inf;

figure("Position",[100,100,300,800]); hold all

data_temp = cell(max_ind,6);
labels = ["B38","I2","B8ab","B6B9","B3","B4B5","Force"];
ylims = [0,20; 0,20; 0,40; 0,50; 0,10; 0,20; 0,300];


for i=1:max_ind
    for j = 1:3
        t_ = data{:,6*(i-1) + 2*(j-1) + 1};
        x_ = data{:,6*(i-1) + 2*(j-1) + 2};
        
        nan_val = isnan(t_);
        t_(nan_val) = [];
        x_(nan_val) = [];
    
        [~,uniqueIdx] = unique(t_);
        t_ = t_(uniqueIdx);
        x_ = x_(uniqueIdx);
        if i==7
            x_([1,end]) = mean(x_([1,end])); % enforcing continuity between cycles
        end
        
        data_temp{i,1+2*(j-1)} = t_;
        data_temp{i,2+2*(j-1)} = x_;

        if min(t_)<t_min
            t_min = min(t_);
        end

        if max(t_)>t_max
            t_max = max(t_);
        end
    end
end

data_interp = cell(max_ind,1);

for i=1:max_ind
    data_ = cell(3,1);
    subplot(7,1,i); hold all
    for j = 1:3
        t_ = data_temp{i,1+2*(j-1)}; t_ = (t_ - t_min)/(t_max - t_min);
        x_ = data_temp{i,2+2*(j-1)}; 
        if i==7
            t_ = [t_; t_+1; t_+2];
            x_ = [x_; x_; x_];
        end

%         data_interp{i,1} = t;
        data_{j,1} = interp1(t_,x_,t,"spline",0);
        
        plot(t,data_{j,1},lines{j})
    end
    data_interp{i} = data_;
    xlim([0,1])
    ylim(ylims(i,:))
    ylabel(labels(i))
end

save(outfile,"data_interp",'-mat')