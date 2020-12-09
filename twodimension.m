clear all
close all

%4 cables, 2 dimensiones (plano), 3 GDL (plano xy, rotacion z)
%sin rozamiento, con efecto de la gravedad.
m = 10; %kg
I_zz = 0.6; 
%% generacion de trayectorias
time=0:0.1:20;

[pos1x,vel1x,ace1x]=gen_coeff(5,6.5,0,0,1.5,-4,0,6);
[pos2x,vel2x,ace2x]=gen_coeff(6.5,3,0,0,-4,2,6,11);
pos2x(:,1)= [];
vel2x(:,1)= [];
ace2x(:,1)= []; %le quito el primer valor para concatenar despues todas las variables en una sola
[pos3x,vel3x,ace3x]=gen_coeff(3,5,0,0,2,0,11,20);
pos3x(:,1)= [];
vel3x(:,1)= [];
ace3x(:,1)= [];

[pos1y,vel1y,ace1y]=gen_coeff(5,2,0,0,-1,5,0,6);
[pos2y,vel2y,ace2y]=gen_coeff(2,6,0,0,5,0.8,6,11);
pos2y(:,1)= [];
vel2y(:,1)= [];
ace2y(:,1)= []; %le quito el primer valor para concatenar despues todas las variables en una sola
[pos3y,vel3y,ace3y]=gen_coeff(6,5,0,0,0.8,0,11,20);
pos3y(:,1)= [];
vel3y(:,1)= [];
ace3y(:,1)= [];

% se puede establecer un limite de -pi/4:pi/4 para la posicion angular a la
% espera del estudio sobre las limitaciones reales
[pos_ang1,vel_ang1,ace_ang1]=gen_coeff(0,0.4,0,0,0.5,-0.2,0,6); %radianes
[pos_ang2,vel_ang2,ace_ang2]=gen_coeff(0.4,-0.2,0,0,-0.2,0.4,6,11);
pos_ang2(:,1)= [];
vel_ang2(:,1)= [];
ace_ang2(:,1)= []; %le quito el primer valor para concatenar despues todas las variables en una sola
[pos_ang3,vel_ang3,ace_ang3]=gen_coeff(-0.2,0,0,0,0.4,0,11,20);
pos_ang3(:,1)= [];
vel_ang3(:,1)= [];
ace_ang3(:,1)= [];

posx=transpose(cat(2,pos1x,pos2x,pos3x));
velx=transpose(cat(2,vel1x,vel2x,vel3x));
acex=transpose(cat(2,ace1x,ace2x,ace3x));

posy=transpose(cat(2,pos1y,pos2y,pos3y));
vely=transpose(cat(2,vel1y,vel2y,vel3y));
acey=transpose(cat(2,ace1y,ace2y,ace3y));

pos_ang=transpose(cat(2,pos_ang1,pos_ang2,pos_ang3));
vel_ang=transpose(cat(2,vel_ang1,vel_ang2,vel_ang3));
ace_ang=transpose(cat(2,ace_ang1,ace_ang2,ace_ang3));

pos(:,1) = posx;
pos(:,2) = posy;
pos(:,3) = pos_ang;

vel(:,1) = velx;
vel(:,2) = vely;
vel(:,3) = vel_ang;

ace(:,1) = acex;
ace(:,2) = acey;
ace(:,3) = ace_ang;

clear('pos1x','pos1y','pos_ang1','pos2x','pos2y','pos_ang2','pos3x','pos3y','pos_ang3','vel1x','vel1y','vel_ang1','vel2x','vel2y','vel_ang2','vel3x','vel3y','vel_ang3','ace1x','ace1y','ace_ang1','ace2x','ace2y','ace_ang2','ace3x','ace3y','ace_ang3','ace_ang','acex','acey','pos_ang','posx','posy','velx','vely','vel_ang')

%% definicion espacio de trabajo
%plano 10x10 m

A = [0 0; 0 10; 10 10; 10 0];

% posición de la plataforma inicial [5 5]
% supongamos una plataforma cuadrada de 1 m de lado
% anclaje cables plataforma, vector desde el centro

B = [0.5 -0.5; 0.5 0.5; -0.5 0.5; -0.5 -0.5];
for i=1:length(time)
Q = [cos(pos(i,3)) -sin(pos(i,3)) ; sin(pos(i,3)) cos(pos(i,3))]; %matriz de rotacion

l = repmat(transpose(pos(i,1:2)),1,4) + Q*transpose(B) - transpose(A); 

L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3)) norm(l(:,4))];

J = zeros(4,3);

J(:,1:2) = [transpose(l)];

J(1,3) = prod_vect(Q*transpose(B(1,:)),transpose(pos(i,1:2))-transpose(A(1,:)));
J(2,3) = prod_vect(Q*transpose(B(2,:)),transpose(pos(i,1:2))-transpose(A(2,:)));
J(3,3) = prod_vect(Q*transpose(B(3,:)),transpose(pos(i,1:2))-transpose(A(3,:)));
J(4,3) = prod_vect(Q*transpose(B(4,:)),transpose(pos(i,1:2))-transpose(A(4,:)));

J(1,:) = J(1,:)/L(1);
J(2,:) = J(2,:)/L(2);
J(3,:) = J(3,:)/L(3);
J(4,:) = J(4,:)/L(4);

W = -transpose(J);

%% analisis dinamico

%M*ace+C*vel+fgravedad=W*t(1,4) caso con fuerzas de coriolis(REVISAR) y gravedad


M = [m 0 0; 0 m 0; 0 0 I_zz];

t_min = 10; %Newtons
t_max = 1000; %Newtons

f_g = [0 -m*9.81 0];

%% prgramacion lineal

    A_des = [-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %matriz de desilguadad Ax<=B
    B_des = [-t_min; -t_min; -t_min; -t_min; t_max; t_max; t_max; t_max];%matriz de desigualdad
    t_opt = zeros(4,length(time)); %inicializo vector de tensiones optimas
    Aeq = W; %matriz de igualdad
    Beq = M*transpose(ace(i,:)); %matriz de igualdad 
    t_opt = linprog([1 1 1 1],A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
    t_final(:,i)= t_opt;
end

plot(time,t_final)
%C = [0 0 (omega x Komega)] nose como deben ponerse estos terminos







