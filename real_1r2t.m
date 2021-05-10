clear all
close all

%% generación de trayectoria deseada

tf = 20;
t0=0;
dt = 0.1;
time=t0:dt:tf;

[posx,velx,acex]=gen_coeff(5,8,0,0,0.5,0,t0,tf,dt);
[posy,vely,acey]=gen_coeff(5,3,0,0,-0.05,0,t0,tf,dt);
[pos_angz,vel_angz,ace_angz]=gen_coeff(0,-0.2,0,0,0-0.03,0,t0,tf,dt);

pos = [posx; posy; pos_angz];
vel = [velx; vely; vel_angz];
ace = [acex; acey; ace_angz];

clear('posx','posy','posz','velx','vely','velz','acex','acey','acez','pos_angx','pos_angy','pos_angz','vel_angx','vel_angy','vel_angz','ace_angx','ace_angy','ace_angz');

%% definicion de parametros invariable y inicializacion

A = [0 0; 0 10; 10 10; 10 0; 5 10]';
B = [-0.5 0.5; -0.5 -0.5; 0.5 -0.5; 0.5 0.5; 0 0.5]';
V = [-0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5]';

t_min = 1;
t_max = 500;

m = 5;
I_zz = 1;
f_g = [0; -m*9.81; 0];

t_linprog = zeros(5,length(time));
t_quadprog = zeros(5,length(time));
t_pinv = zeros(5,length(time));

pos_ode(:,1) = pos(:,1); %cond. ini para los vectores de la trayectoria que seran calculados con ode45
vel_ode(:,1) = vel(:,1);
ace_ode(:,1) = ace(:,1);

%% bucle de calculo - propagación en tiempo

for i=1:length(time)
    
    Q = [cos(pos_ode(3,i)) -sin(pos_ode(3,i)); sin(pos_ode(3,i)) cos(pos_ode(3,i))];
    
    l = repmat(pos_ode(1:2,i),1,5)+Q*B-A;
    
    L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3))...
         norm(l(:,4)) norm(l(:,5))];
    
    J = [(l(:,1:5)./L)' ((l(1,1:5)/L(1,1:5))*(-B(1,1:5)*sin(pos_ode(3,i))-B(2,1:5)*cos(pos_ode(3,i)))...
         +(l(2,1:5)/L)*(B(2,1:5)*cos(pos_ode(3,i))-B(2,1:5)*sin(pos_ode(3,i))))'];
    W = -J';
    
    M = [m*eye(2) zeros(2,1); zeros(1,2) I_zz];
       
   
    %programacion lineal
    A_des = [-eye(5); eye(5)]; %matriz de desilguadad Ax<=B
    B_des = [-t_min*ones(5,1); t_max*ones(5,1)];%matriz de desigualdad
    Aeq = W; %matriz de igualdad
    Beq = M*(ace_ode(:,i)/2)-f_g; %matriz de igualdad EL VALOR DE ACE DEBERIA SER RECALCULADO DESDE EL VALOR REAL DE X(i) hasta la posición ideal X(i+1)
    t_l = linprog(ones(1,5),A_des,B_des,Aeq,Beq);
    if isempty(t_l) == 1
        t_linprog(:,i) = zeros(5,1);
    else
        t_linprog(:,i)= t_l;
    end

    %programacion cuadratica
    A_des = [-eye(5); eye(5)]; %matriz de desilguadad Ax<=B
    B_des = [-t_min*ones(5,1); t_max*ones(5,1)];%matriz de desigualdad
    Aeq = W; %matriz de igualdad
    Beq = M*ace_ode(:,i)-f_g; %matriz de igualdad 
    t_q = quadprog(eye(5),zeros(5,1),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
    if isempty(t_l) == 1
        t_quadprog(:,i) = zeros(5,1);
    else
        t_quadprog(:,i)= t_q;
    end
    
    %inversa de penrose
    t_pinv(:,i) = pinv(W)*(M*ace_ode(:,i) - f_g);
    
    %% ODE45
    % cambiar t_linprog, t_quadprog, t_pinv segun convenga
    ladoderecho = M\(W*t_l+f_g); 
    
    if i<length(time)    %para que no calcule en el instante final
        [t_sol,z] = ode45(@(t,z) odeFunction2(t,z,ladoderecho), [time(i) time(i+1)], [pos_ode(:,i); vel_ode(:,i)]);

        pos_ode(:,i+1) = z(length(z),1:3); %cojo el ultimo valor, que es el que equivale al instante i+1
        vel_ode(:,i+1) = z(length(z),4:6);
        ace_ode(:,i+1) = ladoderecho;
    end

    
    cuadrado = [(pos(1:2,i)+Q*V(:,1))' (pos(1:2,i)+Q*V(:,2))';
                (pos(1:2,i)+Q*V(:,2))' (pos(1:2,i)+Q*V(:,3))';
                (pos(1:2,i)+Q*V(:,3))' (pos(1:2,i)+Q*V(:,4))';
                (pos(1:2,i)+Q*V(:,4))' (pos(1:2,i)+Q*V(:,1))'];


    cables = [(pos(1:2,i)+Q*B(:,1))' A(:,1)';
              (pos(1:2,i)+Q*B(:,2))' A(:,2)';
              (pos(1:2,i)+Q*B(:,3))' A(:,3)';
              (pos(1:2,i)+Q*B(:,4))' A(:,4)';
              (pos(1:2,i)+Q*B(:,5))' A(:,5)'];


    figure(1)
    for k = 1:4
        plot(cuadrado(k,1:2:3),cuadrado(k,2:2:4),'k');
        hold on
    end

    for j = 1:5
        plot(cables(j,1:2:4),cables(j,2:2:4))
        hold on

    end
    axis([0 10 0 10 0 10])
    hold off
    grid on

end