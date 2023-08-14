function fi = Sensor_Inviscid(u,v,w,dx,dy,dz)
% Epsilon
Epsilon = 1E-10;

% Div U
[du_x,du_y,du_z] = CentralDerivative_d1_2ndOrder(u);
[dv_x,dv_y,dv_z] = CentralDerivative_d1_2ndOrder(v);
[dw_x,dw_y,dw_z] = CentralDerivative_d1_2ndOrder(w);

DivU       = du_x./dx + dv_y./dy + dw_z./dz;
DivU2      = DivU.^2;

% Vorticity
Vorticity_x = dw_y./dy - dv_z./dz;
Vorticity_y = du_z./dz - dw_x./dx;
Vorticity_z = dv_x./dx - du_y./dy;

Vort       = sqrt(Vorticity_x.^2 + Vorticity_y.^2 + Vorticity_z.^2);
Vort2      = Vort.^2;

% Sensor activation
fi = DivU2./(DivU2 + Vort2 + Epsilon);


end