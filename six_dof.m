clear all
close all

%supongamos un espacio de trabajo (10x10x10) m


%% definicion de la trayectoria deseada
tf = 20;
t0=0;
time=0:0.1:tf;    
        
[posx,velx,acex]=gen_coeff(5,3,0,0,-0.2,0,t0,tf);
[posy,vely,acey]=gen_coeff(5,8,0,0,0.2,0,t0,tf);
[posz,velz,acez]=gen_coeff(5,6,0,0,0.1,0,t0,tf);
[pos_angx,vel_angx,ace_angx]=gen_coeff(0,0.1,0,0,0.05,0,t0,tf);
[pos_angy,vel_angy,ace_angy]=gen_coeff(0,0.2,0,0,0.04,0,t0,tf);
[pos_angz,vel_angz,ace_angz]=gen_coeff(0,-0.1,0,0,-0.005,0,t0,tf);

pos = [posx; posy; posz; pos_angx; pos_angy; pos_angz];
vel = [velx; vely; velz; vel_angx; vel_angy; vel_angz];
ace = [acex; acey; acez; ace_angx; ace_angy; ace_angz];

clear('posx','posy','posz','velx','vely','velz','acex','acey','acez','pos_angx','pos_angy','pos_angz','vel_angx','vel_angy','vel_angz','ace_angx','ace_angy','ace_angz');

% figure(1)
% plot3(pos(1,:),pos(2,:),pos(3,:));
% axis([0 10 0 10 0 10])
% grid on
% hold on

for i=1:length(time)



%% cinematica y parametros geometricos

Q = [cos(pos(5,i))*cos(pos(6,i)) -cos(pos(5,i))*sin(pos(6,i)) sin(pos(5,i));
cos(pos(4,i))*sin(pos(6,i))+sin(pos(4,i))*sin(pos(5,i))*cos(pos(6,i)) cos(pos(4,i))*cos(pos(6,i))-sin(pos(4,i))*sin(pos(5,i))*sin(pos(6,i)) -sin(pos(4,i))*cos(pos(5,i));
sin(pos(4,i))*sin(pos(6,i))-cos(pos(4,i))*sin(pos(5,i))*cos(pos(6,i)) sin(pos(4,i))*cos(pos(6,i))+cos(pos(4,i))*sin(pos(5,i))*sin(pos(6,i)) cos(pos(4,i))*cos(pos(5,i))];
    

%cables unidos a la plataforma como REELAX8-PC (pag 63) 
%plataforma de 1 m de lado
V = [-0.5 0.5 -0.5; -0.5 -0.5 -0.5; 0.5 -0.5 -0.5 ; 0.5 0.5 -0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 0.5];
A = [0 0 0; 0 10 0; 10 10 0; 10 0 0; 0 0 10; 0 10 10; 10 10 10; 10 0 10];
B = [-0.5 0.5 -0.05; -0.05 -0.5 -0.5; 0.05 -0.5 -0.5 ; 0.5 0.5 -0.05; -0.5 0.5 0.05; -0.05 -0.5 0.5; 0.05 -0.5 0.5; 0.5 0.5 0.05];
l = repmat(pos(1:3,i),1,8)+Q*transpose(B)-transpose(A);

L = [norm(l(:,1:8))];

J = [transpose(l(1:3,:)./L) transpose(cross(Q*transpose(B),repmat(pos(1:3,i),1,8)-transpose(A))./L)];

W = -transpose(J);

%% dinamica
t_min = 40;
t_max = 100000;

m = 10; %kg
I_p = [1.667 0 0; 0 1.667 0; 0 0 1.667];    %se inicializa tensor de inercia de la plataforma

K = Q*I_p*Q';

M = [m*eye(3) zeros(3,3);
    zeros(3,3) K];

f_g = [0; -m*9.81; 0; 0; 0; 0];

C = [zeros(3,1); cross(vel(4:6,i),K*vel(4:6,i))];

% M*ace(:,i)+ C + f_g = W*t 

%programacion lineal
A_des = [-eye(8); eye(8)]; %matriz de desilguadad Ax<=B
B_des = [-t_min*ones(8,1); t_max*ones(8,1)];%matriz de desigualdad
t_opt = zeros(8,length(time)); %inicializo vector de tensiones optimas
Aeq = W; %matriz de igualdad
Beq = M*ace(:,i)+C+f_g; %matriz de igualdad 
t_opt = linprog(ones(1,8),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
t_final(:,i)= t_opt; 

%programacion cuadratica
% A_des = [-eye(8); eye(8)]; %matriz de desilguadad Ax<=B
% B_des = [-t_min*ones(8,1); t_max*ones(8,1)];%matriz de desigualdad
% t_opt = zeros(8,length(time)); %inicializo vector de tensiones optimas
% Aeq = W; %matriz de igualdad
% Beq = M*ace(:,i)+C+f_g; %matriz de igualdad 
% t_opt = quadprog(eye(8),zeros(8,1),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
% t_final_quad(:,i)= t_opt;

%inversa de penrose


%% plotear en tiempo real
% 
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

if i == 2
    W_p = W;
    pos_p = pos(:,2);
end

end

%% validacion. ode45 

% M*ace + C*vel = W*t - f_g

tspan = [0 0.1];
u0 = [pos(:,1)' vel(:,1)'];

[t,x] = ode45(@(t,u) odeFunction(t,u,I_p,Q,M,t_opt,f_g),tspan,u0);




    