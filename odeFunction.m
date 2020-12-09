function du = odeFunction(t,u,I_p,t_opt,W_p,f_g)

x = u(1:6);
v = u(7:12);

Q = [cos(x(5))*cos(x(6)) -cos(x(5))*sin(x(6)) sin(x(5));
cos(x(4))*sin(x(6))+sin(x(4))*sin(x(5))*cos(x(6)) cos(x(4))*cos(x(6))-sin(x(4))*sin(x(5))*sin(x(6)) -sin(x(4))*cos(x(5));
sin(x(4))*sin(x(6))-cos(x(4))*sin(x(5))*cos(x(6)) sin(x(4))*cos(x(6))+cos(x(4))*sin(x(5))*sin(x(6)) cos(x(4))*cos(x(5))];

M = [m*eye(3) zeros(3,3); zeros(3,3) Q*I_p*Q'];


du = [v; inv(M)*(W_p*t_opt - f_g - cross(v(4:6),(Q*I_p*Q')*v(4:6)))];