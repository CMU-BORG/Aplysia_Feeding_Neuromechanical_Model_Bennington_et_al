function Ft = tang_pt_sys(Xt,Xd,Xcom,R)
% system of quadratic equations to find the tangent point on circle that
% passes through Xd. Solve using numerical solver to get tangent point

dorsal = Xd(2)>Xcom(2); % determines if we are dealing with the dorsal or ventral anchor

a = Xd - Xcom; % vector between circle center and anchor point
a_mag = a(1)^2 + a(2)^2;

Ft = [0,0];

% Equation 1: equality of angle from tangent line and line to anchor
LHS1 = (Xd(2) - Xt(2))*(Xt(2) - Xcom(2));
RHS1 = -(Xd(1) - Xt(1))*(Xt(1) - Xcom(1));

Ft(1) = 1*(LHS1 - RHS1);

% Equation 2: length of line between tangent point and anchor point
% LHS2 = Xt(1)^2 + Xt(2)^2 + 2*(a(1) - Xd(1))*Xt(1) + 2*(a(2) - Xd(2))*Xt(2);
% RHS2 = a_mag + R^2 + 2*(a(1)*Xcom(1) + a(2)*Xcom(2)) - (Xd(1)^2 + Xd(2)^2);

LHS2 = (Xt(2)-Xcom(2))^2 + (Xt(1)-Xcom(1))^2;
RHS2 = R^2;

Ft(2) = 1*(LHS2 - RHS2);

% addition of penalty terms to force solution to the back of the
% odontophore
penalty_dorsal =  100*( (Xt(2)<=Xcom(2)) + (Xt(1)>Xd(1)) );
penalty_ventral =  100*( (Xt(2)>=Xcom(2)) + 100*(Xt(1)>Xd(1))  );
penalty = dorsal*penalty_dorsal + (1-dorsal)*penalty_ventral; % + norm(Xt - Xd);
Ft = Ft*(1+penalty);

end
