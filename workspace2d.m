
close all

pos = zeros(5,2601);
k=1;
for i=0:0.2:10
    for j=0:0.2:10
        pos(1,k)=i;
        pos(2,k)=j;
        k=k+1;
    end
end

% %prueba para ver la singularidad
% k=1;
% pos = zeros(5,length(0:0.01:10));
% for i=0:0.01:10
%     pos(2,k)=4;
%     pos(1,k)=i;
%     k=k+1;
% end
% %fin prueba

psi = 0.5;
ace = [0; 0; 0];

A = [0 0; 0 10; 10 10; 10 0; 5 10]';

B = [-0.5 0.5; -0.5 -0.5; 0.5 -0.5; 0.5 0.5; 0 0.5]';

V = [-0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5]';
    
Q = [cos(psi) -sin(psi); sin(psi) cos(psi)];

t_min = 0;
t_max = 500;
m = 50; %kg
f_g = [0; -9.8*m; 0];
I_zz = 8.33; %kg*m^2

t_penrose = zeros(5,length(pos));
t_linprog = zeros(5,length(pos));
t_quadprog = zeros(5,length(pos));

M = [m*eye(2) zeros(2,1); zeros(1,2) I_zz];

singu = zeros(5,length(pos));

counter = 0;
%k=332;
for k=1:length(pos)

    l = repmat(pos(1:2,k),1,5)+Q*B-A;

    L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3)) norm(l(:,4)) norm(l(:,5))];

    J = [transpose(l(:,1:5)./L) transpose((l(1,1:5)/L(1,1:5))*(-B(1,1:5)*sin(psi)-B(2,1:5)*cos(psi))...
        +(l(2,1:5)/L)*(B(2,1:5)*cos(psi)-B(2,1:5)*sin(psi)))];

    W = -transpose(J);
    
    %recoleccion datos W(3,1) para comprobar singularidad
    
    singu(:,k)=W(3,:);
    
    %fin prueba recoleccion
    
    t_pinv = lsqminnorm(W,M*ace-f_g); %es lo mismo que pseudoinversa de penrose
    t_penrose(:,k) = t_pinv; 
   
    A_des = [-eye(5); eye(5)]; %matriz de desilguadad Ax<=B
    B_des = [-t_min*ones(5,1); t_max*ones(5,1)];%matriz de desigualdad
    Aeq = W; %matriz de igualdad
    Beq = M*ace-f_g; %matriz de igualdad                                     
    
    % t(5x1) = t_pinv(5x1) + N(5x2)*lambda(2x1)
    % tmin-t_pinv < N * lambda < tmin-t_pinv
   
    t_l = linprog(ones(1,5),A_des,B_des,Aeq,Beq);
    if isempty(t_l) == 1
        t_linprog(:,k) = zeros(5,1);
        pos(3,k)=0;
    else
        t_linprog(:,k)= t_l;
        pos(3,k)=1;
    end

    
    t_q = quadprog(eye(5),zeros(5,1),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
    if isempty(t_q) == 1
        t_quadprog(:,k) = zeros(5,1);
        pos(4,k)=0;
    else
        t_quadprog(:,k)= t_q;
        pos(4,k)=1;
    end

    
    alpha_c(1) = pi+atan(abs(l(2,1))/abs(l(1,1)));
    alpha_c(2) = (pi/2)+atan(abs(l(1,2))/abs(l(2,2)));
    alpha_c(3) = atan(abs(l(2,3))/abs(l(1,3)));
    alpha_c(4) = pi*(3/2)+atan(abs(l(1,4))/abs(l(2,4)));
    
    alpha_c_grados = alpha_c*(180/pi);
    
    alpha_p(1) = psi+(3/2)*pi;
    alpha_p(2) = psi+0.5*pi;
    alpha_p(3) = psi+0.5*pi;
    alpha_p(4) = psi+(3/2)*pi;
    
    angulo_plataforma_grados = alpha_p*(180/pi);
    
    if alpha_p(1)<=alpha_c(1) || alpha_p(3)<=alpha_c(3) || alpha_c(2)<=alpha_p(2) || alpha_c(4)<=alpha_p(4)
        pos(5,k)=0;
    else
        pos(5,k)=1;
    end
    
    if pos(5,k) == 1 && pos(4,k) == 1
        counter = counter+1;
    end
                           
end

porcentaje = counter/length(pos);
%% bucle de interferencias

figure(1)
grid on
for k=1:length(pos)
    if pos(5,k) == 0
        plot(pos(1,k),pos(2,k),'b*');
        hold on
    elseif pos(3,k) == 0 && pos(4,k) == 0 
        plot(pos(1,k),pos(2,k),'r.','MarkerSize',10);
        hold on
    elseif pos(3,k) == 0 && pos(4,k) == 1
        plot(pos(1,k),pos(2,k),'go','MarkerSize',4);
        hold on    
    else
        plot(pos(1,k),pos(2,k),'g.','MarkerSize',10);
    end
end
xlabel('x')
ylabel('y')

% figure(2)
% plot(pos(1,:),singu(1,:),'r.-')
% xlabel('x')
% ylabel('W(3,1)')
