close all

N = null(W);

coord_x = -2000:1:2000;

v = -2000:1:2000;  
[x, y] = meshgrid(v);  
conditions = (t_min-t_pinv(1)<x*N(1,1)+y*N(1,2)) &...
             (t_min-t_pinv(2)<x*N(2,1)+y*N(2,2)) &...
             (t_min-t_pinv(3)<x*N(3,1)+y*N(3,2)) &...
             (t_min-t_pinv(4)<x*N(4,1)+y*N(4,2)) &...
             (t_min-t_pinv(5)<x*N(5,1)+y*N(5,2)) &...
             (t_min-t_pinv(6)<x*N(6,1)+y*N(6,2)) &...
             (t_min-t_pinv(7)<x*N(7,1)+y*N(7,2)) &...
             (t_min-t_pinv(8)<x*N(8,1)+y*N(8,2)) &...
             (t_max-t_pinv(1)>x*N(1,1)+y*N(1,2)) &...
             (t_max-t_pinv(2)>x*N(2,1)+y*N(2,2)) &...
             (t_max-t_pinv(3)>x*N(3,1)+y*N(3,2)) &...
             (t_max-t_pinv(4)>x*N(4,1)+y*N(4,2)) &...
             (t_max-t_pinv(5)>x*N(5,1)+y*N(5,2)) &...
             (t_max-t_pinv(6)>x*N(6,1)+y*N(6,2)) &...
             (t_max-t_pinv(7)>x*N(7,1)+y*N(7,2)) &...
             (t_max-t_pinv(8)>x*N(8,1)+y*N(8,2));
cond = zeros(length(v)); % Initialize
cond(conditions) = NaN;
surf(x, y, cond)
view(0,90)
t_bari = t_pinv+N*[200; 100];
hold on

y_1min=(t_min-t_pinv(1)-coord_x*N(1,1))/N(1,2);
y_2min=(t_min-t_pinv(2)-coord_x*N(2,1))/N(2,2);
y_3min=(t_min-t_pinv(3)-coord_x*N(3,1))/N(3,2);
y_4min=(t_min-t_pinv(4)-coord_x*N(4,1))/N(4,2);
y_5min=(t_min-t_pinv(5)-coord_x*N(5,1))/N(5,2);
y_6min=(t_min-t_pinv(6)-coord_x*N(6,1))/N(6,2);
y_7min=(t_min-t_pinv(7)-coord_x*N(7,1))/N(7,2);
y_8min=(t_min-t_pinv(8)-coord_x*N(8,1))/N(8,2);
y_1max=(t_max-t_pinv(1)-coord_x*N(1,1))/N(1,2);
y_2max=(t_max-t_pinv(2)-coord_x*N(2,1))/N(2,2);
y_3max=(t_max-t_pinv(3)-coord_x*N(3,1))/N(3,2);
y_4max=(t_max-t_pinv(4)-coord_x*N(4,1))/N(4,2);
y_5max=(t_max-t_pinv(5)-coord_x*N(5,1))/N(5,2);
y_6max=(t_max-t_pinv(6)-coord_x*N(6,1))/N(6,2);
y_7max=(t_max-t_pinv(7)-coord_x*N(7,1))/N(7,2);
y_8max=(t_max-t_pinv(8)-coord_x*N(8,1))/N(8,2);

plot(coord_x,y_1min,'r--')
plot(coord_x,y_1max,'r--')
plot(coord_x,y_2min,'r--')
plot(coord_x,y_2max,'r--')
plot(coord_x,y_3min,'r--')
plot(coord_x,y_3max,'r--')
plot(coord_x,y_4min,'r--')
plot(coord_x,y_4max,'r--')
plot(coord_x,y_5min,'r--')
plot(coord_x,y_5max,'r--')
plot(coord_x,y_6min,'r--')
plot(coord_x,y_6max,'r--')
plot(coord_x,y_7min,'r--')
plot(coord_x,y_7max,'r--')
plot(coord_x,y_8min,'r--')
plot(coord_x,y_8max,'r--')
axis([-2000 2000 -2000 2000])

xlabel('\lambda_1')
ylabel('\lambda_2')
ace_calculada = M\(W*t_bari+f_g);