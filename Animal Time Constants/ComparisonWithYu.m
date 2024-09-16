%% Fitting Double First Order Filter Activation Function to Yu et al 1999 activation function
% Last Updated: MB 6/10/2024

% I2 parameters
tau = 2.45;     % activation time constant
beta = 0.703;   % relaxation time constant correction

% preallocating memory
tau_fit = zeros(9,1);

figure("Position",[100,100,940,640],"Color","w"); hold all

% for 9 different stimulation durations
for i=1:9
    subplot(3,3,i); hold all
    ta = 1;     % time to turn on stimulation
    tb = 3*i;   % time to turn off stimulation
    tend = 45;  % duration of simulation
    
    u = @(t) (t>ta).*(t<=tb);   % stimulation function
    
    da_prime_dt = @(t,a) (1/tau)*( u(t) - (beta + (1-beta)*u(t))*a );   % activation nonlinear governing equation
    
    % solving for level of activation
    [t,a] = ode23(da_prime_dt,linspace(0,tend,1000),0,odeset("RelTol",1e-5,"AbsTol",1e-4));
    
    % error between the nonlinear activation solution and double first
    % order filter solution
    err = @(x) sum( ((a - Tension(t,x(1),x(1)/x(2),ta,tb))).^2 );
    taus = fminsearch(err,[tau,beta]); % optimize parameters, starting from the parameters in Yu 1999

    % saving params
    tau_fit(i) = taus(1);
    beta_fit(i) = taus(2);

    % plotting
    fill([ta,ta,tb,tb],[0,1,1,0],0.6*[1,1,1],"EdgeColor","none","FaceAlpha",0.2,"HandleVisibility","off")

    plot(t,a,'-r',"LineWidth",1,"DisplayName","Yu 1999 Model")
    plot(t,Tension(t,taus(1),taus(1)/taus(2),ta,tb),'-b',"LineWidth",1,"DisplayName","Double-First-Order Filter Model")
    if i==2
        legend("NumColumns",2)
    end
    xlim([0,tend])
    xlabel("Time [s]")
    ylabel("Activation [ ]")
    set(gca,"FontSize",12)
end
saveas(gcf,"Figures\ComparisonWithYu.svg",'svg')

figure("Position",[100,100,940,300],"color","w");
subplot(1,2,1); hold all
plot(3*(1:9) - 1,tau_fit,'ok','MarkerFaceColor','k',"MarkerSize",4,"DisplayName","Individual Fits")
yline(median(tau_fit),'--b',"Displayname","Median","LineWidth",1)
legend("NumColumns",2)
xlabel("Duration of Stimulation [s]")
ylabel("Optimized \tau_{on}")
ylim([0,2])
set(gca,"FontSize",15)

subplot(1,2,2); hold all
plot(3*(1:9) - 1,beta_fit,'ok','MarkerFaceColor','k',"MarkerSize",4,"DisplayName","Individual Fits")
yline(median(beta_fit),'--g',"Displayname","Median","LineWidth",1)
legend("NumColumns",2)
xlabel("Duration of Stimulation [s]")
ylabel("Optimized \tau_{on}/\tau_{off}")
ylim([0,2])
set(gca,"FontSize",15)
saveas(gcf,"Figures\ConvergenceWithYu.svg",'svg')


fprintf("Double First Order Filter Time Constants: tau = %.3f, beta = %.3f\n", median(tau_fit), median(beta_fit))


function T = Tension(t,tau_on,tau_off,t_on,t_off)


T_rise = @(t) ( 1 - exp(-(t-t_on)/tau_on) .* (1 + (t-t_on)/tau_on) );
T_peak = T_rise(t_off);
A_peak = 1 - exp(-(t_off-t_on)/tau_on);

T_fall = @(t) exp(-(t-t_off)/tau_off) .* (T_peak + A_peak*(t-t_off)/tau_off);

T = (t>t_on).*(t<=t_off).*T_rise(t) + (t>t_off).*T_fall(t);

end