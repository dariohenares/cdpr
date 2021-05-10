N = null(W);

v = -2000:1:2000;  % plotting range from -5 to 5
[x, y] = meshgrid(v);  % get 2-D mesh for x and y
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

t_bari = t_pinv+N*[600; 0];