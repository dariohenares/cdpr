close all
clear all
%% generacion de trayectorias
tf = 10;
t0=0;
dt = 0.1;
time=t0:dt:tf;    

pos = zeros(6,length(time));
vel = zeros(6,length(time));
ace = zeros(6,length(time));

waypoint_t = [5 5 5; 3 5 7; 4 5 6; 3.5 3.5 3.5; 5 5 5]';
waypoint_r = [0 0 0; 0 0 0; 0 0.25 -0.1; 0.2 -0.1 0.25; 0 0 0]';
times = [t0 2.5 5 7.5 tf];
[pos(1:3,:), vel(1:3,:), ace(1:3,:)] = quinticpolytraj(waypoint_t,times,time);
[pos(4:6,:), vel(4:6,:), ace(4:6,:)] = quinticpolytraj(waypoint_r,times,time);

%% definicion de parametros invariables e inicializacion

V = [-0.5 -0.5 -0.5; -0.5 0.5 -0.5; 0.5 0.5 -0.5 ; 0.5 -0.5 -0.5; -0.5 -0.5 0.5; -0.5 0.5 0.5; 0.5 0.5 0.5 ; 0.5 -0.5 0.5]';
A = [0 0 0; 0 10 0; 10 10 0; 10 0 0; 0 0 10; 0 10 10; 10 10 10; 10 0 10]';
B = [-0.5 -0.05 0.5; -0.5 0.05 0.5; 0.5 0.05 0.5; 0.5 -0.05 0.5; -0.05 -0.5 -0.5; -0.05 0.5 -0.5; 0.05 0.5 -0.5; 0.05 -0.5 -0.5]'; %nuevo


%B = [-0.5 0.5 -0.05; -0.05 -0.5 -0.5; 0.05 -0.5 -0.5 ; 0.5 0.5 -0.05; -0.5 0.5 0.05; -0.05 -0.5 0.5; 0.05 -0.5 0.5; 0.5 0.5 0.05]; %antiguo
%B = [-0.5 -0.5 -0.5; -0.5 0.5 -0.5; 0.5 0.5 -0.5 ; 0.5 -0.5 -0.5; -0.5 -0.5 0.5; -0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 -0.5 0.5];
%B = [-0.5 -0.05 -0.5; -0.5 0.05 -0.5; 0.5 0.05 -0.5 ; 0.5 -0.05 -0.5; -0.05 -0.5 0.5; -0.05 0.5 0.5; 0.05 0.5 0.5; 0.05 -0.5 0.5];


t_min = 0;
t_max = 2000;
lb = ones(8,1)*t_min;
ub = ones(8,1)*t_max; %valores limite para problema de optimizacion


m = 100; %kg
I_p = [16.67 0 0; 0 16.67 0; 0 0 16.67]; %tensor inercia
f_g = [0; 0; -m*9.81; 0; 0; 0]; %fuerza gravedad
t_linprog = zeros(8,length(time));
t_quadprog = zeros(8,length(time));
t_penrose = zeros(8,length(time));

pos_ode = zeros(6,length(time)); 
vel_ode = zeros(6,length(time)); 
ace_ode = zeros(6,length(time));

pos_ode(:,1) = pos(:,1); %cond. ini para los vectores de la trayectoria que seran calculados con ode45
vel_ode(:,1) = vel(:,1);
ace_ode(:,1) = ace(:,1);

dist = zeros(4,length(time));
dist_p = zeros(16,length(time));

cubo_plot = [];
cables_plot = [];
colors = [rgb('Red'); rgb('Goldenrod'); rgb('Green'); rgb('Blue'); rgb('Purple'); rgb('DimGrey'); rgb('SaddleBrown'); rgb('Aqua')];

%% Blucle de cálculo - progagación en tiempo
tic
for i=1:length(time)

    Q = [cos(pos(5,i))*cos(pos(6,i)) -cos(pos(5,i))*sin(pos(6,i)) sin(pos(5,i));
        cos(pos(4,i))*sin(pos(6,i))+sin(pos(4,i))*sin(pos(5,i))*cos(pos(6,i))...
        cos(pos(4,i))*cos(pos(6,i))-sin(pos(4,i))*sin(pos(5,i))*sin(pos(6,i))...
        -sin(pos(4,i))*cos(pos(5,i));
        sin(pos(4,i))*sin(pos(6,i))-cos(pos(4,i))*sin(pos(5,i))*cos(pos(6,i))...
        sin(pos(4,i))*cos(pos(6,i))+cos(pos(4,i))*sin(pos(5,i))*sin(pos(6,i))...
        cos(pos(4,i))*cos(pos(5,i))]; %matriz de rotacion
    
    
    l = repmat(pos(1:3,i),1,8)+Q*B-A; %matriz de longitudes cables (xyz) 
    
    L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3)) norm(l(:,4))...
         norm(l(:,5)) norm(l(:,6)) norm(l(:,7)) norm(l(:,8))]; %vector modulo longitudes cables

    J = [(l(1:3,:)./L)' (cross(Q*B,repmat(pos(1:3,i),1,8)-A)./L)']; %matriz jacobiana
    W = -J'; %matriz estructura
    K = Q*I_p*Q'; %tensor de inercia rotado
    M = [m*eye(3) zeros(3,3); zeros(3,3) K]; %matriz M de la ec. dinamica
    C = [zeros(3,1); cross(vel(4:6,i),K*vel(4:6,i))]; %vector c de ec. dinamica
    
    %% interferencias
    %interferencias entre pares de cables
    
    dist(1,i) = line_to_line_distance(A(:,1),l(:,1)/L(1),A(:,5),-l(:,5)/L(5));
    dist(2,i) = line_to_line_distance(A(:,2),-l(:,2)/L(2),A(:,6),-l(:,6)/L(6));
    dist(3,i) = line_to_line_distance(A(:,3),l(:,3)/L(3),A(:,7),-l(:,7)/L(7));
    dist(4,i) = line_to_line_distance(A(:,4),-l(:,4)/L(4),A(:,8),-l(:,8)/L(8));
    
    %interferencias entre cable y aristas
    dist_p(1,i) = line_to_line_distance(A(:,1),-l(:,1)/L(1),pos(1:3,i)+Q*V(:,1),Q*[0;1;0]);
    dist_p(2,i) = line_to_line_distance(A(:,1),l(:,1)/L(1),pos(1:3,i)+Q*V(:,1),Q*[0;0;1]);
    dist_p(3,i) = line_to_line_distance(A(:,2),-l(:,2)/L(2),pos(1:3,i)+Q*V(:,2),Q*[0;1;0]);
    dist_p(4,i) = line_to_line_distance(A(:,2),-l(:,2)/L(2),pos(1:3,i)+Q*V(:,2),Q*[0;0;1]);
    dist_p(5,i) = line_to_line_distance(A(:,3),l(:,3)/L(3),pos(1:3,i)+Q*V(:,3),Q*[0;1;0]);
    dist_p(6,i) = line_to_line_distance(A(:,3),l(:,3)/L(3),pos(1:3,i)+Q*V(:,3),Q*[0;0;1]);
    dist_p(7,i) = line_to_line_distance(A(:,4),l(:,4)/L(4),pos(1:3,i)+Q*V(:,4),Q*[0;1;0]);
    dist_p(8,i) = line_to_line_distance(A(:,4),-l(:,4)/L(4),pos(1:3,i)+Q*V(:,4),Q*[0;0;1]);
    dist_p(9,i) = line_to_line_distance(A(:,5),-l(:,5)/L(5),pos(1:3,i)+Q*V(:,5),Q*[1;0;0]);
    dist_p(10,i) = line_to_line_distance(A(:,5),-l(:,5)/L(5),pos(1:3,i)+Q*V(:,5),Q*[0;0;1]);
    dist_p(11,i) = line_to_line_distance(A(:,6),l(:,6)/L(6),pos(1:3,i)+Q*V(:,6),Q*[1;0;0]);
    dist_p(12,i) = line_to_line_distance(A(:,6),l(:,6)/L(6),pos(1:3,i)+Q*V(:,6),Q*[0;0;1]);
    dist_p(13,i) = line_to_line_distance(A(:,7),l(:,7)/L(7),pos(1:3,i)+Q*V(:,7),Q*[1;0;0]);
    dist_p(14,i) = line_to_line_distance(A(:,7),-l(:,7)/L(7),pos(1:3,i)+Q*V(:,7),Q*[0;0;1]);
    dist_p(15,i) = line_to_line_distance(A(:,8),-l(:,8)/L(8),pos(1:3,i)+Q*V(:,8),Q*[1;0;0]);
    dist_p(16,i) = line_to_line_distance(A(:,8),l(:,8)/L(8),pos(1:3,i)+Q*V(:,8),Q*[0;0;1]);


    %% optimizacion lineal y cuadratica
    Aeq = W; 
    Beq = M*ace(:,i)+C-f_g; %ecuacion de equilibrio de fuerzas
   
    t_l = linprog(ones(1,8),[],[],Aeq,Beq,lb,ub);
    if isempty(t_l) == 1
        t_linprog(:,i) = zeros(8,1);
    else
        t_linprog(:,i)= t_l;
    end 
 
    t_q = quadprog(eye(8),zeros(8,1),[],[],Aeq,Beq,lb,ub); %t_opt  vector de tensiones de los cables
    if isempty(t_q) == 1
        t_quadprog(:,i) = zeros(8,1);
    else
        t_quadprog(:,i)= t_q; 
    end

    t_pinv = pinv(W)*(M*ace(:,i) + C - f_g);
    t_penrose(:,i) = t_pinv;
    N = null(W);
    
%     lambda_l = N\(t_l-t_pinv);
%     lambda_q = N\(t_q-t_pinv);

    
    %% almacenamiento de variables para plotear

    cubo = [(pos(1:3,i)+Q*V(:,1))' (pos(1:3,i)+Q*V(:,2))';
            (pos(1:3,i)+Q*V(:,2))' (pos(1:3,i)+Q*V(:,3))';
            (pos(1:3,i)+Q*V(:,3))' (pos(1:3,i)+Q*V(:,4))';
            (pos(1:3,i)+Q*V(:,4))' (pos(1:3,i)+Q*V(:,1))';
            (pos(1:3,i)+Q*V(:,5))' (pos(1:3,i)+Q*V(:,6))';
            (pos(1:3,i)+Q*V(:,6))' (pos(1:3,i)+Q*V(:,7))';
            (pos(1:3,i)+Q*V(:,7))' (pos(1:3,i)+Q*V(:,8))';
            (pos(1:3,i)+Q*V(:,8))' (pos(1:3,i)+Q*V(:,5))';
            (pos(1:3,i)+Q*V(:,1))' (pos(1:3,i)+Q*V(:,5))';
            (pos(1:3,i)+Q*V(:,2))' (pos(1:3,i)+Q*V(:,6))';
            (pos(1:3,i)+Q*V(:,3))' (pos(1:3,i)+Q*V(:,7))';
            (pos(1:3,i)+Q*V(:,4))' (pos(1:3,i)+Q*V(:,8))'];
          
    cubo_plot = cat(2,cubo_plot,cubo);
        
    cables = [(pos(1:3,i)+Q*B(:,1:8)); A(:,1:8)]';
    
    cables_plot = cat(2,cables_plot,cables);
    
end
toc
%% bucle 2 ODE
for i=1:length(time)
    
    Q_ode = [cos(pos_ode(5,i))*cos(pos_ode(6,i)) -cos(pos_ode(5,i))*sin(pos_ode(6,i)) sin(pos_ode(5,i));
    cos(pos_ode(4,i))*sin(pos_ode(6,i))+sin(pos_ode(4,i))*sin(pos_ode(5,i))*cos(pos_ode(6,i))...
    cos(pos_ode(4,i))*cos(pos_ode(6,i))-sin(pos_ode(4,i))*sin(pos_ode(5,i))*sin(pos_ode(6,i))...
    -sin(pos_ode(4,i))*cos(pos_ode(5,i));
    sin(pos_ode(4,i))*sin(pos_ode(6,i))-cos(pos_ode(4,i))*sin(pos_ode(5,i))*cos(pos_ode(6,i))...
    sin(pos_ode(4,i))*cos(pos_ode(6,i))+cos(pos_ode(4,i))*sin(pos_ode(5,i))*sin(pos_ode(6,i))...
    cos(pos_ode(4,i))*cos(pos_ode(5,i))]; %matriz de rotacion
        
    l_ode = repmat(pos_ode(1:3,i),1,8)+Q_ode*B-A; %matriz de longitudes cables (xyz) 
    
    L_ode = [norm(l_ode(:,1)) norm(l_ode(:,2)) norm(l_ode(:,3)) norm(l_ode(:,4))...
         norm(l_ode(:,5)) norm(l_ode(:,6)) norm(l_ode(:,7)) norm(l_ode(:,8))]; %vector modulo longitudes cables

    J_ode = [(l_ode(1:3,:)./L_ode)' (cross(Q_ode*B,repmat(pos_ode(1:3,i),1,8)-A)./L_ode)']; %matriz jacobiana
    W_ode = -J_ode'; %matriz estructura
    K_ode = Q_ode*I_p*Q_ode'; %tensor de inercia rotado
    M_ode = [m*eye(3) zeros(3,3); zeros(3,3) K_ode]; %matriz M de la ec. dinamica
    C_ode = [zeros(3,1); cross(vel_ode(4:6,i),K*vel_ode(4:6,i))]; %vector c de ec. dinamica
    
    ladoderecho = M_ode\(W_ode*t_quadprog(:,i)+f_g-C_ode); 
    
    if i<length(time)    %para que no calcule en el instante final
        [t_sol,z] = ode45(@(t,z) odeFunction(t,z,ladoderecho), [time(i) time(i+1)], [pos_ode(:,i); vel_ode(:,i)]);

        pos_ode(:,i+1) = z(length(z),1:6); %cojo el ultimo valor, que es el que equivale al instante i+1
        vel_ode(:,i+1) = z(length(z),7:12);
        ace_ode(:,i+1) = ladoderecho;
    end
end 



 %% PLOTS
 
% figure(1)
% set(gcf, 'Position',  [100, 100, 700, 600])
% subplot(2,1,1)
% for i=1:8
%     plot(time(:), t_penrose(i,:),'Color',colors(i,:),'LineWidth',0.8)
%     hold on
% end
% ylabel('Tensiones (N)')
% axis([0 20 -500 700])
% legend('1','2','3','4','5','6','7','8')
% subplot(2,1,2)
% for i=1:8
%     plot(time(:), t_quadprog(i,:),'Color',colors(i,:),'LineWidth',0.8)
%     hold on
% end
% ylabel('Tensiones (N)')
% xlabel('Tiempo (s)')
% axis([0 20 0 1200])
% legend('1','2','3','4','5','6','7','8')
% subplot(2,1,2)

%% PLOTS
myvideo = VideoWriter('Prueba_defensa');
myvideo.FrameRate = 10;
open(myvideo)

cont = 0;
figure(1)
set(gcf, 'Position',  [20, 50, 1500, 650])
for i=1:6:6*length(time)
    subplot(1,2,1)
    for k = 1:12
        plot3(cubo_plot(k,i:3:i+3),cubo_plot(k,i+1:3:i+4),cubo_plot(k,i+2:3:i+5),'k');
        hold on
    end

    for j = 1:8
        plot3(cables_plot(j,i:3:i+3),cables_plot(j,i+1:3:i+4),cables_plot(j,i+2:3:i+5),'Color',colors(j,:),'LineWidth',0.8)
        hold on
    end
   
    axis([0 10 0 10 0 10])
    hold off
    grid off
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view (-65,18)
    
    %base fija
    line([0 0],[0 10],[0 0],'color','k')
    line([0 0],[0 0],[0 10],'color','k')
    line([0 0],[10 10],[0 10],'color','k')
    line([0 0],[0 10],[10 10],'color','k')
    line([10 10],[0 10],[0 0],'color','k')
    line([10 10],[0 0],[0 10],'color','k')
    line([10 10],[10 10],[0 10],'color','k')
    line([10 10],[0 10],[10 10],'color','k')
    line([0 10],[10 10],[10 10],'color','k')
    line([0 10],[0 0],[10 10],'color','k')
    line([0 10],[0 0],[0 0],'color','k')
    line([0 10],[10 10],[0 0],'color','k')
    
    cont=cont+1;
    
    txt1=['Distancia minima entre cable 1 y 5: ' num2str(dist(1,cont)) ' m '];
    txt2=['Distancia minima entre cable 2 y 6: ' num2str(dist(2,cont)) ' m '];
    txt3=['Distancia minima entre cable 3 y 7: ' num2str(dist(3,cont)) ' m '];
    txt4=['Distancia minima entre cable 4 y 8: ' num2str(dist(4,cont)) ' m '];
    text1_p=['cable 1 - arista opuesta: ',num2str(dist_p(1,cont)),' m ; cable 1 - arista lateral: ',num2str(dist_p(2,cont)),' m'];
    text2_p=['cable 2 - arista opuesta: ',num2str(dist_p(3,cont)),' m ; cable 2 - arista lateral: ',num2str(dist_p(4,cont)),' m'];
    text3_p=['cable 3 - arista opuesta: ',num2str(dist_p(5,cont)),' m ; cable 3 - arista lateral: ',num2str(dist_p(6,cont)),' m'];
    text4_p=['cable 4 - arista opuesta: ',num2str(dist_p(7,cont)),' m ; cable 4 - arista lateral: ',num2str(dist_p(8,cont)),' m'];
    text5_p=['cable 5 - arista opuesta: ',num2str(dist_p(9,cont)),' m ; cable 5 - arista lateral: ',num2str(dist_p(10,cont)),' m'];
    text6_p=['cable 6 - arista opuesta: ',num2str(dist_p(11,cont)),' m ; cable 6 - arista lateral: ',num2str(dist_p(12,cont)),' m'];
    text7_p=['cable 7 - arista opuesta: ',num2str(dist_p(13,cont)),' m ; cable 7 - arista lateral: ',num2str(dist_p(14,cont)),' m'];
    text8_p=['cable 8 - arista opuesta: ',num2str(dist_p(15,cont)),' m ; cable 8 - arista lateral: ',num2str(dist_p(16,cont)),' m'];
    hold off
    
    subplot(1,2,2)
    for j=1:8
        plot(time(1:cont),t_quadprog(j,1:cont),'Color',colors(j,:),'LineWidth',1)
        hold on
    end
    text(1,1600,{txt1 txt2 txt3 txt4 text1_p text2_p text3_p text4_p text5_p text6_p text7_p text8_p},'FontSize',8)
    axis([0 10 0 2000])
    grid on
    xlabel('Tiempo (s)')
    ylabel('Tensiones (N)')
    legend('1','2','3','4','5','6','7','8')
    frame = getframe(gcf);
    writeVideo(myvideo, frame);
    hold off
end


figure(2)
set(gcf, 'Position',  [100, 100, 1400, 600])
tit = ["x","y","z","\phi","\theta","\psi"];
ylab = ["Pos (m)","Pos (m)","Pos (m)","Pos (rad)","Pos (rad)","Pos (rad)"];
for i=1:6
    subplot(2,3,i)
    plot(time,pos(i,:)); hold on; plot(time,pos_ode(i,:));
    grid on
    xlabel('Tiempo (s)')
    ylabel(ylab(i))
    title(tit(i))
    legend('Teorica','Real')
end

figure(3)
set(gcf, 'Position',  [100, 100, 1400, 600])
tit = ["v_x","v_y","v_z","\omega_x","\omega_y","\omega_z"];
ylab = ["Vel (m/s)","Vel (m/s)","Vel (m/s)","Vel (rad/s)","Vel (rad/s)","Vel (rad/s)"];
for i=1:6
    subplot(2,3,i)
    plot(time,vel(i,:)); hold on; plot(time,vel_ode(i,:));
    grid on
    xlabel('Tiempo (s)')
    ylabel(ylab(i))
    title(tit(i))
    legend('Teorica','Real')
end

figure(4)
set(gcf, 'Position',  [100, 100, 1400, 600])
tit = ["a_x","a_y","a_z","\alpha_x","\alpha_y","\alpha_z"];
ylab = ["Ace (m/s^2)","Ace (m/s^2)","Ace (m/s^2)","Ace (rad/s^2)","Ace (rad/s^2)","Ace (rad/s^2)",];
for i=1:6
    subplot(2,3,i)
    plot(time,ace(i,:)); hold on; plot(time,ace_ode(i,:));
    grid on
    xlabel('Tiempo (s)')
    ylabel(ylab(i))
    title(tit(i))
    legend('Teorica','Real')
end


figure(5)
error_pos = abs(pos(1:6,:)-pos_ode(1:6,:));
plot(time,error_pos);
grid on
legend('x','y','z','\phi','\theta','\psi')
title('Error de posición acumulado')
xlabel('Tiempo (s)')
ylabel('Error de posición (m/rad)')

% 
% figure(6)
% plot3(pos(1,:),pos(2,:),pos(3,:),'r.');
% axis([0 10 0 10 0 10])
% grid on
% hold on


