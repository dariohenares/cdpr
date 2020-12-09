clear all
close all
m = 20; %kg
%dos cables
%suponemos que el centro de gravedad de la plataforma concuerda con el
%centro geométrico, ademas las cuerdas no generan ningun momento ya que
%la prolongación de las cuerdas pasan por CG.

%caso: tensiones como variable de salida (no rozamiento)

J = [1; -1];
W = -(transpose(J));
t = zeros(2,1);

%% generacion de trayectorias

%polinomio de 5º orden: x(t) = a*t^5+b*t^4+c*t^3+d*t^2+e*t+f
time=0:0.1:20;

[pos1,vel1,ace1]=gen_coeff(0,5,0,0,3,0,0,6);
[pos2,vel2,ace2]=gen_coeff(5,-8,0,0,-2,1,6,11);
pos2(:,1)= [];
vel2(:,1)= [];
ace2(:,1)= []; %le quito el primer valor para concatenar despues todas las variables en una sola
[pos3,vel3,ace3]=gen_coeff(-8,0,0,0,1,0,11,20);
pos3(:,1)= [];
vel3(:,1)= [];
ace3(:,1)= [];


pos=transpose(cat(2,pos1,pos2,pos3));
vel=transpose(cat(2,vel1,vel2,vel3));
ace=transpose(cat(2,ace1,ace2,ace3));

figure(1)
subplot(3,1,1)
plot(time,pos)
ylabel('Posición (m)')
subplot(3,1,2)
plot(time,vel)
ylabel('Velocidad (m/s)')
subplot(3,1,3)
plot(time,ace)
ylabel('Aceleracion (m/s^2)')

%% analisis dinámico: m*ace = W*t

%condiciones t
%t>0, tmin<t<tmax

t_min = 10; %Newtons
t_max = 500; %Newtons

%m*[ace]201x1 = -[t1]201x1 + [t2]201x1

%% programacion lineal (SIMPLEX)    
    A = [-1 0; 0 -1; 1 0; 0 1]; %matriz de desilguadad Ax<=B
    B = [-t_min; -t_min; t_max; t_max];%matriz de desigualdad
    t_opt = zeros(2,length(time)); %inicializo vector de tensiones optimas
for i = 1:length(time)
    Aeq = W; %matriz de igualdad
    Beq = ace(i)*m; %matriz de igualdad -t1+t2=ace*m
    t_opt(:,i) = linprog(W,A,B,Aeq,Beq); %t_opt  vector de tensiones de los cables
end

figure(2)
plot(time,t_opt(1,:),'r','DisplayName','t_1')
hold on;
plot(time,t_opt(2,:),'b','DisplayName','t_2')
ylabel('N')
legend
axis([0 20 0 100])










