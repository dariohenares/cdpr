clear all
close all

%% definicion de la trayectoria deseada
tf = 20;
t0=0;
dt = 0.1;
time=t0:dt:tf;    
        
[posx,velx,acex]=gen_coeff(5,3,0,0,-0.2,0,t0,tf,dt);
[posy,vely,acey]=gen_coeff(5,8,0,0,0.2,0,t0,tf,dt);
[posz,velz,acez]=gen_coeff(5,6,0,0,0.1,0,t0,tf,dt);
[pos_angx,vel_angx,ace_angx]=gen_coeff(0,0.1,0,0,0.05,0,t0,tf,dt);
[pos_angy,vel_angy,ace_angy]=gen_coeff(0,0.2,0,0,0.04,0,t0,tf,dt);
[pos_angz,vel_angz,ace_angz]=gen_coeff(0,-0.1,0,0,-0.005,0,t0,tf,dt);

pos = [posx; posy; posz; pos_angx; pos_angy; pos_angz];
vel = [velx; vely; velz; vel_angx; vel_angy; vel_angz];
ace = [acex; acey; acez; ace_angx; ace_angy; ace_angz];

clear('posx','posy','posz','velx','vely','velz','acex','acey','acez','pos_angx','pos_angy','pos_angz','vel_angx','vel_angy','vel_angz','ace_angx','ace_angy','ace_angz');

% figure(1)
% plot3(pos(1,:),pos(2,:),pos(3,:));
% axis([0 10 0 10 0 10])
% grid on
% hold on

%% definicion de parametros invariables e inicializacion

V = [-0.5 0.5 -0.5; -0.5 -0.5 -0.5; 0.5 -0.5 -0.5 ; 0.5 0.5 -0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 0.5];
A = [0 0 0; 0 10 0; 10 10 0; 10 0 0; 0 0 10; 0 10 10; 10 10 10; 10 0 10];
B = [-0.5 0.5 -0.05; -0.05 -0.5 -0.5; 0.05 -0.5 -0.5 ; 0.5 0.5 -0.05; -0.5 0.5 0.05; -0.05 -0.5 0.5; 0.05 -0.5 0.5; 0.5 0.5 0.05];

t_min = 1;
t_max = 500;
m = 10; %kg
I_p = [1.667 0 0; 0 1.667 0; 0 0 1.667]; %tensor inercia
f_g = [0; -m*9.81; 0; 0; 0; 0]; %fuerza gravedad
t_linprog = zeros(8,length(time));
t_quadprog = zeros(8,length(time));

pos_ode = zeros(6,length(time)); 
vel_ode = zeros(6,length(time)); 
ace_ode = zeros(6,length(time));

pos_ode(:,1) = pos(:,1); %cond. ini para los vectores de la trayectoria que seran calculados con ode45
vel_ode(:,1) = vel(:,1);
ace_ode(:,1) = ace(:,1);

%% Blucle de cálculo - progagación en tiempo

for i=1:length(time)

    Q = [cos(pos(5,i))*cos(pos(6,i)) -cos(pos(5,i))*sin(pos(6,i)) sin(pos(5,i));
        cos(pos(4,i))*sin(pos(6,i))+sin(pos(4,i))*sin(pos(5,i))*cos(pos(6,i))...
        cos(pos(4,i))*cos(pos(6,i))-sin(pos(4,i))*sin(pos(5,i))*sin(pos(6,i))...
        -sin(pos(4,i))*cos(pos(5,i));
        sin(pos(4,i))*sin(pos(6,i))-cos(pos(4,i))*sin(pos(5,i))*cos(pos(6,i))...
        sin(pos(4,i))*cos(pos(6,i))+cos(pos(4,i))*sin(pos(5,i))*sin(pos(6,i))...
        cos(pos(4,i))*cos(pos(5,i))]; %matriz de rotacion
    
    
    l = repmat(pos(1:3,i),1,8)+Q*B'-A'; %matriz de longitudes cables (xyz) 
    
    L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3)) norm(l(:,4))...
         norm(l(:,5)) norm(l(:,6)) norm(l(:,7)) norm(l(:,8))]; %vector modulo longitudes cables

    J = [(l(1:3,:)./L)' (cross(Q*B',repmat(pos(1:3,i),1,8)-A')./L)']; %matriz jacobiana

    W = -J'; %matriz estructura

    K = Q*I_p*Q'; %tensor de inercia rotado

    M = [m*eye(3) zeros(3,3); zeros(3,3) K]; %matriz M de la ec. dinamica

    C = [zeros(3,1); cross(vel(4:6,i),K*vel(4:6,i))]; %vector c de ec. dinamica

    % M*ace(:,i)+ C = W*t + f_g

    %programacion lineal
    A_des = [-eye(8); eye(8)]; %matriz de desilguadad Ax<=B
    B_des = [-t_min*ones(8,1); t_max*ones(8,1)];%matriz de desigualdad
    Aeq = W; %matriz de igualdad
    Beq = M*ace(:,i)+C-f_g; %matriz de igualdad 
    t_l = linprog(ones(1,8),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
    if isempty(t_l) == 1
        t_linprog(:,i) = zeros(8,1);
    else
        t_linprog(:,i)= t_l;
    end 

    %programacion cuadratica
    A_des = [-eye(8); eye(8)]; %matriz de desilguadad Ax<=B
    B_des = [-t_min*ones(8,1); t_max*ones(8,1)];%matriz de desigualdad
    Aeq = W; %matriz de igualdad
    Beq = M*ace(:,i)+C-f_g; %matriz de igualdad 
    t_q = quadprog(eye(8),zeros(8,1),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
    if isempty(t_l) == 1
        t_quadprog(:,i) = zeros(8,1);
    else
        t_quadprog(:,i)= t_q; 
    end

    %inversa de penrose
    t_pinv = pinv(W)*(M*ace(:,i) + C - f_g);


    % plotear en tiempo real

    % cubo = [transpose(pos(1:3,i)+Q*transpose(V(1,:))) transpose(pos(1:3,i)+Q*transpose(V(2,:)));
    %         transpose(pos(1:3,i)+Q*transpose(V(2,:))) transpose(pos(1:3,i)+Q*transpose(V(3,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(3,:))) transpose(pos(1:3,i)+Q*transpose(V(4,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(4,:))) transpose(pos(1:3,i)+Q*transpose(V(1,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(1,:))) transpose(pos(1:3,i)+Q*transpose(V(5,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(5,:))) transpose(pos(1:3,i)+Q*transpose(V(6,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(6,:))) transpose(pos(1:3,i)+Q*transpose(V(7,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(7,:))) transpose(pos(1:3,i)+Q*transpose(V(8,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(8,:))) transpose(pos(1:3,i)+Q*transpose(V(5,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(2,:))) transpose(pos(1:3,i)+Q*transpose(V(6,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(3,:))) transpose(pos(1:3,i)+Q*transpose(V(7,:)))
    %         transpose(pos(1:3,i)+Q*transpose(V(4,:))) transpose(pos(1:3,i)+Q*transpose(V(8,:)))];
    % 
    % cables = [transpose(pos(1:3,i)+Q*transpose(B(1:8,:))) A(1:8,:)];
    % 
    % 
    % figure(1)
    % set(gcf, 'Position',  [100, 100, 1400, 600])
    % subplot(1,2,1)
    % view([-37.5 -30])
    % for k = 1:12
    %     plot3(cubo(k,1:3:4),cubo(k,2:3:5),cubo(k,3:3:6),'k');
    %     hold on
    % end
    % 
    % for j = 1:8
    %     plot3(cables(j,1:3:4),cables(j,2:3:5),cables(j,3:3:6))
    %     hold on
    % 
    % end
    % axis([0 10 0 10 0 10])
    % hold off
    % grid on
    % 
    % subplot(1,2,2)
    % plot(time(1:i),t_final(:,1:i))
    % axis([0 20 0 1000])
    % hold off
    % grid on
    
    %% ODE45
    % cambiar t_linprog, t_quadprog, t_pinv segun convenga
    ladoderecho = M\(W*t_l+f_g-C); 
    
    if i<length(time)    %para que no calcule en el instante final
        [t_sol,z] = ode45(@(t,z) odeFunction(t,z,ladoderecho), [time(i) time(i+1)], [pos_ode(:,i); vel_ode(:,i)]);

        pos_ode(:,i+1) = z(length(z),1:6); %cojo el ultimo valor, que es el que equivale al instante i+1
        vel_ode(:,i+1) = z(length(z),7:12);
        ace_ode(:,i+1) = ladoderecho;
    end

end

%% PLOTS

figure(1)
set(gcf, 'Position',  [100, 100, 1400, 600])
tit = ["x","y","z","\phi","\theta","\psi"];
ylab = ["Pos (m)","Pos (m)","Pos (m)","Pos (rad)","Pos (rad)","Pos (rad)"];
for i=1:6
    subplot(2,3,i)
    plot(time,pos(i,:)); hold on; plot(time,pos_ode(i,:));
    grid on
    xlabel('Time (s)')
    ylabel(ylab(i))
    title(tit(i))
end

figure(2)
set(gcf, 'Position',  [100, 100, 1400, 600])
tit = ["v_x","v_y","v_z","\omega_x","\omega_y","\omega_z"];
ylab = ["Vel (m/s)","Vel (m/s)","Vel (m/s)","Vel (rad/s)","Vel (rad/s)","Vel (rad/s)"];
for i=1:6
    subplot(2,3,i)
    plot(time,vel(i,:)); hold on; plot(time,vel_ode(i,:));
    grid on
    xlabel('Time (s)')
    ylabel(ylab(i))
    title(tit(i))
    legend('Teorica','Real')
end

figure(3)
set(gcf, 'Position',  [100, 100, 1400, 600])
tit = ["a_x","a_y","a_z","\alpha_x","\alpha_y","\alpha_z"];
ylab = ["Ace (m/s^2)","Ace (m/s^2)","Ace (m/s^2)","Ace (rad/s^2)","Ace (rad/s^2)","Ace (rad/s^2)",];
for i=1:6
    subplot(2,3,i)
    plot(time,ace(i,:)); hold on; plot(time,ace_ode(i,:));
    grid on
    xlabel('Time (s)')
    ylabel(ylab(i))
    title(tit(i))
    legend('Teorica','Real')
end


figure(4)
error_pos = abs(pos(1:6,:)-pos_ode(1:6,:));
plot(time,error_pos);
grid on
legend('x','y','z','\phi','\theta','\psi')
title('Accumulated positional error')
xlabel('Time (s)')
ylabel('Position error (m)')
