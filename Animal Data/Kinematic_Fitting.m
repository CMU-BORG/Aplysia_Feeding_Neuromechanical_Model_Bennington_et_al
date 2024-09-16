function [mean_func,std_func] = Kinematic_Fitting(mean_data,mean_p1std_data)

datasets = {mean_data,mean_p1std_data};
interps = cell(2,1);

for i=1:2
    ds = datasets{i};
    nans = isnan(ds(:,1));
    ds(nans,:) = [];
       
    x = ds(:,1);
    x = [x;x+x(end);x+2*x(end)];
    y = ds(:,2);
    y = [y;y;y];
    
    [xData, yData] = prepareCurveData( x, y );
    
    % Set up fittype and options.
    ft = fittype( 'fourier8' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [mean(y) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2*pi];
    
    % Fit model to data.
    [fitresult, ~] = fit( xData, yData, ft, opts );
    interps{i} = @(t) fitresult(t);
end

mean_func = interps{1};
std_func = @(t) interps{2}(t) - interps{1}(t);

end