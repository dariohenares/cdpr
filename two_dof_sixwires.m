close all

%% generación de trayectoria deseada

tf = 20;
t0=0;
dt = 0.1;
time=t0:dt:tf;

pos = zeros(3,length(time));
vel = zeros(3,length(time));
ace = zeros(3,length(time));

r = 3.5;
vert = zeros(5,2);
center = [5 4];

for i=1:5
    vert(i,1)=r*cos(i*2*pi/5);
    vert(i,2)=r*sin(i*2*pi/5);
end

rot = -pi/10;
Q = [cos(rot) -sin(rot); sin(rot) cos(rot)];
vert(:,1:2)=vert(:,1:2)*Q;
vert(:,1:2) = vert(:,1:2)+center(1:2);

waypoint = [vert(1,1) vert(1,2); vert(3,1) vert(3,2); vert(5,1) vert(5,2); vert(2,1) vert(2,2); vert(4,1) vert(4,2); vert(1,1) vert(1,2)]';
times = [t0 2 7 12 14 tf];
[pos(1:2,:), vel(1:2,:), ace(1:2,:)] = quinticpolytraj(waypoint,times,time);

[pos(3,:), vel(3,:), ace(3,:)] = quinticpolytraj([0 -0.4 0.1 -0.15 0.4 0],times,time);

%% definicion de parametros invariables y inicializacion

%A = [0 0; 0 10; 10 10; 10 0; 5 10]';
A = [0 0; 0 10; 10 10; 10 0; 0 10; 10 10]';
%B = [-0.5 0.5; -0.5 -0.5; 0.5 -0.5; 0.5 0.5; 0 0.5]';
B = [-0.5 0.5; -0.5 -0.5; 0.5 -0.5; 0.5 0.5; 0 0.5; 0 0.5]';
V = [-0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5]';

t_min = 0;
t_max = 500;
lb = ones(6,1)*t_min;
ub = ones(6,1)*t_max;

m = 50;
I_zz = 1;
f_g = [0; -m*9.81; 0];

t_linprog = zeros(6,length(time));
t_quadprog = zeros(6,length(time));
t_penrose = zeros(6,length(time));

pos_ode(:,1) = pos(:,1); %cond. ini para los vectores de la trayectoria que seran calculados con ode45
vel_ode(:,1) = vel(:,1);
ace_ode(:,1) = ace(:,1);

colors = [rgb('Red'); rgb('Gold'); rgb('Green'); rgb('Blue'); rgb('Purple'); rgb('DimGrey')];

cuadrado_plot = [];
cables_plot = [];
anglebetween = zeros(4,length(time));

%% bucle de calculo - propagación en tiempo

for i=1:length(time)
    
    Q = [cos(pos(3,i)) -sin(pos(3,i)); sin(pos(3,i)) cos(pos(3,i))];
    
    l = repmat(pos(1:2,i),1,6)+Q*B-A;
    
    L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3))...
         norm(l(:,4)) norm(l(:,5)) norm(l(:,6))];
    J = [(l(:,1:6)./L)' ((l(1,1:6)/L(1,1:6))*(-B(1,1:6)*sin(pos(3,i))-B(2,1:6)*cos(pos(3,i)))...
         +(l(2,1:6)/L)*(B(2,1:6)*cos(pos(3,i))-B(2,1:6)*sin(pos(3,i))))'];
    W = -J';
    M = [m*eye(2) zeros(2,1); zeros(1,2) I_zz];
    
    %% interferencias
    
    alpha_c(1) = pi+atan(abs(l(2,1))/abs(l(1,1)));
    alpha_c(2) = (pi/2)+atan(abs(l(1,2))/abs(l(2,2)));
    alpha_c(3) = atan(abs(l(2,3))/abs(l(1,3)));
    alpha_c(4) = pi*(3/2)+atan(abs(l(1,4))/abs(l(2,4)));
    
    alpha_c = alpha_c*(180/pi);
    
    alpha_p(1) = pos(3,i)+(3/2)*pi;
    alpha_p(2) = pos(3,i)+0.5*pi;
    alpha_p(3) = pos(3,i)+0.5*pi;
    alpha_p(4) = pos(3,i)+(3/2)*pi;
    
    alpha_p = alpha_p*(180/pi);
    
    anglebetween(:,i) = [alpha_p(1)-alpha_c(1); alpha_c(2)-alpha_p(2); alpha_p(3)-alpha_c(3); alpha_c(4)-alpha_p(4)];
    
    %% optimizacion
    Aeq = W; 
    Beq = M*ace(:,i)-f_g;  
    
    t_l = linprog(ones(1,6),[],[],Aeq,Beq,lb,ub);
    if isempty(t_l) == 1
        t_linprog(:,i) = zeros(6,1);
    else
        t_linprog(:,i)= t_l;
    end

    t_q = quadprog(eye(6),zeros(6,1),[],[],Aeq,Beq,lb,ub); 
    if isempty(t_l) == 1
        t_quadprog(:,i) = zeros(6,1);
    else
        t_quadprog(:,i)= t_q;
    end
    
    t_pinv = pinv(W)*(M*ace(:,i) - f_g);
    t_penrose(:,i) = t_pinv;
    
    %% ODE45
    ladoderecho = M\(W*t_q+f_g); 
    
    if i<length(time)
        [t_sol,z] = ode45(@(t,z) odeFunction2(t,z,ladoderecho), [time(i) time(i+1)], [pos_ode(:,i); vel_ode(:,i)]);

        pos_ode(:,i+1) = z(length(z),1:3); %cojo el ultimo valor, que es el que equivale al instante i+1
        vel_ode(:,i+1) = z(length(z),4:6);
        ace_ode(:,i+1) = ladoderecho;
    end

      
    %% almacenamiento de variables para plotear
    
    cuadrado = [(pos(1:2,i)+Q*V(:,1))' (pos(1:2,i)+Q*V(:,2))';
                (pos(1:2,i)+Q*V(:,2))' (pos(1:2,i)+Q*V(:,3))';
                (pos(1:2,i)+Q*V(:,3))' (pos(1:2,i)+Q*V(:,4))';
                (pos(1:2,i)+Q*V(:,4))' (pos(1:2,i)+Q*V(:,1))'];
            
    cables = [(pos(1:2,i)+Q*B(:,1))' A(:,1)';
              (pos(1:2,i)+Q*B(:,2))' A(:,2)';
              (pos(1:2,i)+Q*B(:,3))' A(:,3)';
              (pos(1:2,i)+Q*B(:,4))' A(:,4)';
              (pos(1:2,i)+Q*B(:,5))' A(:,5)'
              (pos(1:2,i)+Q*B(:,6))' A(:,6)'];

    cuadrado_plot = cat(2,cuadrado_plot,cuadrado);
    cables_plot = cat(2,cables_plot,cables);
    
    
end

%% Plots
cont = 0;
figure(1)
set(gcf, 'Position',  [100, 100, 1400, 600])
for i=1:4:4*length(time)
      
    subplot(1,2,1)
    for k = 1:4
        plot(cuadrado_plot(k,i:2:i+2),cuadrado_plot(k,i+1:2:i+3),'k');
        hold on
    end

    for j = 1:6
        plot(cables_plot(j,i:2:i+2),cables_plot(j,i+1:2:i+3),'Color',colors(j,:),'LineWidth',0.8)
        hold on
    end
    axis([0 10 0 10 0 10])
    hold off
    grid on
    xlabel('x')
    ylabel('y')
    
    cont=cont+1;
    
    txt1=['Angulo entre cable 1 y plataforma ' num2str(anglebetween(1,cont)) 'º'];
    txt2=['Angulo entre cable 2 y plataforma ' num2str(anglebetween(2,cont)) 'º'];
    txt3=['Angulo entre cable 3 y plataforma ' num2str(anglebetween(3,cont)) 'º'];
    txt4=['Angulo entre cable 4 y plataforma ' num2str(anglebetween(4,cont)) 'º'];
    
    subplot(1,2,2)
    for j=1:6
        plot(time(1:cont),t_quadprog(j,1:cont),'Color',colors(j,:),'LineWidth',0.8);
        hold on
    end
    text(1,650,{txt1 txt2 txt3 txt4},'FontSize',8)
    axis([0 20 0 700])
    xlabel('Tiempo (s)')
    ylabel('Tensiones (N)')
    legend('1','2','3','4','5')
    hold off
    pause(0.05)
end
% 
% figure(2)
% set(gcf, 'Position',  [100, 100, 1400, 600])
% tit = ["x","y","\psi"];
% ylab = ["Pos (m)","Pos (m)","Pos (rad)"];
% for i=1:3
%     subplot(1,3,i)
%     plot(time,pos(i,:)); hold on; plot(time,pos_ode(i,:));
%     grid on
%     xlabel('Time (s)')
%     ylabel(ylab(i))
%     title(tit(i))
% end
% 
% figure(3)
% set(gcf, 'Position',  [100, 100, 1400, 600])
% tit = ["v_x","v_y","\omega_z"];
% ylab = ["Vel (m/s)","Vel (m/s)","Vel (rad/s)"];
% for i=1:3
%     subplot(1,3,i)
%     plot(time,vel(i,:)); hold on; plot(time,vel_ode(i,:));
%     grid on
%     xlabel('Time (s)')
%     ylabel(ylab(i))
%     title(tit(i))
%     legend('Teorica','Real')
% end
% 
% figure(4)
% set(gcf, 'Position',  [100, 100, 1400, 600])
% tit = ["a_x","a_y","\alpha_z"];
% ylab = ["Ace (m/s^2)","Ace (m/s^2)","Ace (rad/s^2)",];
% for i=1:3
%     subplot(1,3,i)
%     plot(time,ace(i,:)); hold on; plot(time,ace_ode(i,:));
%     grid on
%     xlabel('Time (s)')
%     ylabel(ylab(i))
%     title(tit(i))
%     legend('Teorica','Real')
% end
% 
% 
% figure(5)
% error_pos = abs(pos(1:3,:)-pos_ode(1:3,:));
% plot(time,error_pos);
% grid on
% legend('x','y','\psi')
% title('Accumulated positional error')
% xlabel('Time (s)')
% ylabel('Position error (m)')