function [pos,vel,ace]=gen_coeff(pos0,posf,vel0,velf,ace0,acef,t0,tf,dt)
t = t0:dt:tf;

B = [pos0; posf; vel0; velf; ace0; acef];
A = [t0^5 t0^4 t0^3 t0^2 t0 1;
    tf^5 tf^4 tf^3 tf^2 tf 1;
    5*t0^4 4*t0^3 3*t0^2 2*t0 1 0;
    5*tf^4 4*tf^3 3*tf^2 2*tf 1 0;
    20*t0^3 12*t0^2 6*t0 2 0 0;
    20*tf^3 12*tf^2 6*tf 2 0 0];

coeff = A\B;



pos = coeff(1)*t.^5+coeff(2)*t.^4+coeff(3)*t.^3+coeff(4)*t.^2+coeff(5)*t+coeff(6);
vel = 5*coeff(1)*t.^4+4*coeff(2)*t.^3+3*coeff(3)*t.^2+2*coeff(4)*t+coeff(5);
ace = 20*coeff(1)*t.^3+12*coeff(2)*t.^2+6*coeff(3)*t+2*coeff(4);



%x_0=a*t_0^5+b*t_0^4+c*t_0^3+d*t_0^2+e*t_0+f
%x_f=a*t_f^5+b*t_f^4+c*t_f^3+d*t_f^2+e*t_f+f
%x_dot_0=5*a*t_0^4+4*b*t_0^3+3*c*t_0^2+2*d*t_0+e
%x_dot_f=5*a*t_f^4+4*b*t_f^3+3*c*t_f^2+2*d*t_f+e
%x_dot_dot_0=20*a*t_0^3+12*b*t_0^2+6*c*t_0+2*d
%x_dot_dot_f=20*a*t_0^3+12*b*t_0^2+6*c*t_0+2*d


