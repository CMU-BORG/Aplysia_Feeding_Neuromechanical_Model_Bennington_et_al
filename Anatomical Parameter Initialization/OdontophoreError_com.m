function err = OdontophoreError_com(X_data,R1,R2,L,xg,yg,thetag,V0,params)
    dir = split(cd,"\");
    dir = strjoin(dir(1:end-1),"\");
    addpath(dir)

    N = 100;
    THETA = linspace(0,2*pi,N);
    params.phi = GetTangent(R1,R2,L);

    f_vol = @(alpha) V0 - OdontophoreVolume([R1,R2,L,alpha],params);
    alpha = fsolve(f_vol,1,optimset('Display','off'));

    [x_pts,y_pts] = DrawOdontophoreMS([xg,yg,thetag,R1,R2,L,alpha],params);
    [~,solved] = GetTangent(R1,R2,L);

    theta_model = atan2(y_pts - yg,x_pts - xg)';
    R_model = sqrt((x_pts - xg).^2 + (y_pts - yg).^2)';
    theta_model = [theta_model(1:end-1)-2*pi, theta_model, theta_model(2:end)+2*pi];
    [theta_model,u_idx] = unique(theta_model);
    R_model = [R_model(1:end-1), R_model, R_model(2:end)];
    R_model = R_model(u_idx);
    R_model_func = @(x) interp1(theta_model,R_model,x,'spline');
    X_model_func = @(x) [xg;yg] + R_model_func(x).*[cos(x);sin(x)];

    err = sum(sum((X_data(THETA) - X_model_func(THETA)).^2)) + ...
          (1e10)*(1-solved) + ...
          0*(1/R1)^2 + ...
          0*(1/R2)^2 + ...
          100*((2*max(R1,R2) - (R1+R2+L)) > 0.001);
end
