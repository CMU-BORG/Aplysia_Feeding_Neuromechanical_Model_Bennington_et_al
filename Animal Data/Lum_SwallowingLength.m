%{
Length of Swallowed Seaweed
Reproduced from Lum et al 2005

%}
file = "Lum2005_SwallowLength_Raw.csv";
folder = cd + "\Datasets\";

t = linspace(0,1,100);

data_swallow = readmatrix(folder+file,'NumHeaderLines',2);
swallow1 = data_swallow(:,[1,2]);
swallow1(isnan(swallow1(:,1)),:) = [];
swallow1(:,1) = (swallow1(:,1) - min(swallow1(:,1)))/(max(swallow1(:,1)) - min(swallow1(:,1)) );
swallow1(:,2) = (swallow1(:,2) - (swallow1(1,2)));%/(max(swallow1(:,2)) - (swallow1(1,2)) );
swallow1_ = @(t) interp1(swallow1(:,1),swallow1(:,2),t,'spline');

swallow2 = data_swallow(:,[3,4]);
swallow2(isnan(swallow2(:,1)),:) = [];
swallow2(:,1) = (swallow2(:,1) - min(swallow2(:,1)))/(max(swallow2(:,1)) - min(swallow2(:,1)) );
swallow2(:,2) = (swallow2(:,2) - (swallow2(1,2)));%/(max(swallow2(:,2)) - (swallow2(1,2)) );
swallow2_ = @(t) interp1(swallow2(:,1),swallow2(:,2),t,'spline');

swallow3 = data_swallow(:,[5,6]);
swallow3(isnan(swallow3(:,1)),:) = [];
swallow3(:,1) = (swallow3(:,1) - min(swallow3(:,1)))/(max(swallow3(:,1)) - min(swallow3(:,1)) );
swallow3(:,2) = (swallow3(:,2) - (swallow3(1,2)));%/(max(swallow3(:,2)) - (swallow3(1,2)) );
swallow3_ = @(t) interp1(swallow3(:,1),swallow3(:,2),t,'spline');

swallow_mean_data = mean([swallow1_(t);swallow2_(t);swallow3_(t)],1);
swallow_std_data = std([swallow1_(t);swallow2_(t);swallow3_(t)],[],1);
length_swallow_mean = @(t_) interp1(t,swallow_mean_data,t_,'spline');
length_swallow_std = @(t_) interp1(t,swallow_std_data,t_,'spline');

figure; hold all
plot(t,swallow1_(t),'.r','MarkerSize',2);
plot(t,swallow2_(t),'.b','MarkerSize',2);
plot(t,swallow3_(t),'.m','MarkerSize',2);
plot(t,length_swallow_mean(t),'-k','LineWidth',1);
plot(t,length_swallow_mean(t) + length_swallow_std(t),'--k','LineWidth',0.75)
plot(t,length_swallow_mean(t) - length_swallow_std(t),'--k','LineWidth',0.75)
