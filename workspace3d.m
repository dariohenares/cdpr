clear all
close all

pos = zeros(5,17576);
z=1;
for i=0:0.4:10
    for j=0:0.4:10
        for k=0:0.4:10
            pos(1,z)=i;
            pos(2,z)=j;
            pos(3,z)=k;
            z=z+1;
        end 
    end
end

V = [-0.5 -0.5 -0.5; -0.5 0.5 -0.5; 0.5 0.5 -0.5 ; 0.5 -0.5 -0.5; -0.5 -0.5 0.5; -0.5 0.5 0.5; 0.5 0.5 0.5 ; 0.5 -0.5 0.5]';
A = [0 0 0; 0 10 0; 10 10 0; 10 0 0; 0 0 10; 0 10 10; 10 10 10; 10 0 10]';
B = [-0.5 -0.05 0.5; -0.5 0.05 0.5; 0.5 0.05 0.5; 0.5 -0.05 0.5; -0.05 -0.5 -0.5; -0.05 0.5 -0.5; 0.05 0.5 -0.5; 0.05 -0.5 -0.5]';
o = [0.2617993878 0.2617993878 0.2617993878]; %orientatacion dada a la plataforma
%.5235987756
t_min = 0;
t_max = 1000;
m = 100; %kg
f_g = [0; 0; -m*9.81; 0; 0; 0]; %fuerza gravedad
t_linprog = zeros(8,length(pos));
t_quadprog = zeros(8,length(pos));
t_penrose = zeros(8,length(pos));

dist = zeros(1,4);

ace = [0; 0; 0; 0; 0; 0];
I_p = [16.67 0 0; 0 16.67 0; 0 0 16.67];

counter = 0;

%% bucle
%k=3088;
for k=1:length(pos)

    Q = [cos(o(2))*cos(o(3)) -cos(o(2))*sin(o(3)) sin(o(2));
        cos(o(1))*sin(o(3))+sin(o(1))*sin(o(2))*cos(o(3))...
        cos(o(1))*cos(o(3))-sin(o(1))*sin(o(2))*sin(o(3))...
        -sin(o(1))*cos(o(2));
        sin(o(1))*sin(o(3))-cos(o(1))*sin(o(2))*cos(o(3))...
        sin(o(1))*cos(o(3))+cos(o(1))*sin(o(2))*sin(o(3))...
        cos(o(1))*cos(o(2))]; %matriz de rotacion
    
    
    l = repmat(pos(1:3,k),1,8)+Q*B-A; %matriz de longitudes cables (xyz) 
    
    L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3)) norm(l(:,4))...
         norm(l(:,5)) norm(l(:,6)) norm(l(:,7)) norm(l(:,8))]; %vector modulo longitudes cables

    J = [(l(1:3,:)./L)' (cross(Q*B,repmat(pos(1:3,k),1,8)-A)./L)']; %matriz jacobiana

    W = -J'; %matriz estructura
    
    %interferencias cable-cable
    dist(1) = line_to_line_distance(A(:,1),l(:,1)/L(1),A(:,5),-l(:,5)/L(5));
    dist(2) = line_to_line_distance(A(:,2),-l(:,2)/L(2),A(:,6),-l(:,6)/L(6));
    dist(3) = line_to_line_distance(A(:,3),l(:,3)/L(3),A(:,7),-l(:,7)/L(7));
    dist(4) = line_to_line_distance(A(:,4),-l(:,4)/L(4),A(:,8),-l(:,8)/L(8));
    
    %interferencias cable-plataforma
    dist_p(1) = line_to_line_distance(A(:,1),-l(:,1)/L(1),pos(1:3,k)+Q*V(:,1),Q*[0;1;0]);
    dist_p(2) = line_to_line_distance(A(:,1),l(:,1)/L(1),pos(1:3,k)+Q*V(:,1),Q*[0;0;1]);
    dist_p(3) = line_to_line_distance(A(:,2),-l(:,2)/L(2),pos(1:3,k)+Q*V(:,2),Q*[0;1;0]);
    dist_p(4) = line_to_line_distance(A(:,2),-l(:,2)/L(2),pos(1:3,k)+Q*V(:,2),Q*[0;0;1]);
    dist_p(5) = line_to_line_distance(A(:,3),l(:,3)/L(3),pos(1:3,k)+Q*V(:,3),Q*[0;1;0]);
    dist_p(6) = line_to_line_distance(A(:,3),l(:,3)/L(3),pos(1:3,k)+Q*V(:,3),Q*[0;0;1]);
    dist_p(7) = line_to_line_distance(A(:,4),l(:,4)/L(4),pos(1:3,k)+Q*V(:,4),Q*[0;1;0]);
    dist_p(8) = line_to_line_distance(A(:,4),-l(:,4)/L(4),pos(1:3,k)+Q*V(:,4),Q*[0;0;1]);
    dist_p(9) = line_to_line_distance(A(:,5),-l(:,5)/L(5),pos(1:3,k)+Q*V(:,5),Q*[1;0;0]);
    dist_p(10) = line_to_line_distance(A(:,5),-l(:,5)/L(5),pos(1:3,k)+Q*V(:,5),Q*[0;0;1]);
    dist_p(11) = line_to_line_distance(A(:,6),l(:,6)/L(6),pos(1:3,k)+Q*V(:,6),Q*[1;0;0]);
    dist_p(12) = line_to_line_distance(A(:,6),l(:,6)/L(6),pos(1:3,k)+Q*V(:,6),Q*[0;0;1]);
    dist_p(13) = line_to_line_distance(A(:,7),l(:,7)/L(7),pos(1:3,k)+Q*V(:,7),Q*[1;0;0]);
    dist_p(14) = line_to_line_distance(A(:,7),-l(:,7)/L(7),pos(1:3,k)+Q*V(:,7),Q*[0;0;1]);
    dist_p(15) = line_to_line_distance(A(:,8),-l(:,8)/L(8),pos(1:3,k)+Q*V(:,8),Q*[1;0;0]);
    dist_p(16) = line_to_line_distance(A(:,8),l(:,8)/L(8),pos(1:3,k)+Q*V(:,8),Q*[0;0;1]);

    if dist(1)<0 || dist(2)<0 || dist(3)<0 || dist(4)<0 || dist_p(1)<0 || dist_p(2)<0 || dist_p(3)<0 || dist_p(4)<0 ||...
            dist_p(5)<0 || dist_p(6)<0 || dist_p(7)<0 || dist_p(8)<0 || dist_p(9)<0 || dist_p(10)<0 ||...
            dist_p(11)<0 || dist_p(12)<0 || dist_p(13)<0 || dist_p(14)<0 || dist_p(15)<0 || dist_p(16)<0
        pos(5,k)=0;
    else
        pos(5,k)=1;
    end
    
    K = Q*I_p*Q'; %tensor de inercia rotado
    
    M = [m*eye(3) zeros(3,3); zeros(3,3) K]; %matriz M de la ec. dinamica
    
    A_des = [-eye(8); eye(8)]; %matriz de desilguadad Ax<=B
    B_des = [-t_min*ones(8,1); t_max*ones(8,1)];%matriz de desigualdad
    Aeq = W; %matriz de igualdad
    Beq = M*ace-f_g; %matriz de igualdad
    
    % se comenta linprog para que vaya mas rapido el calculo
%     t_l = linprog(ones(1,8),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
%     if isempty(t_l) == 1
%         t_linprog(:,k) = zeros(8,1);
%         pos(4,k) = 0;
%     else
%         t_linprog(:,k)= t_l;
%         pos(4,k) = 1;
%     end 

  
    t_q = quadprog(eye(8),zeros(8,1),A_des,B_des,Aeq,Beq); %t_opt  vector de tensiones de los cables
    if isempty(t_q) == 1
        t_quadprog(:,k) = zeros(8,1);
        pos(4,k) = 0;
    else
        t_quadprog(:,k) = t_q;
        pos(4,k) = 1;
    end

    if pos(5,k) == 1 && pos(4,k) == 1
        counter = counter+1;
    end
    
    
end

porcentaje = counter/length(pos);

figure(1)
for k=1:length(pos)
    if pos(5,k) == 1 && pos(4,k) == 1
        plot3(pos(1,k),pos(2,k),pos(3,k),'g.','MarkerSize',5)
        hold on
    end
end
axis([0 10 0 10 0 10])
xlabel('x')
ylabel('y')
zlabel('z')
grid on

    
figure(2)
for k=1:length(pos)
    if pos(5,k) == 0
        plot3(pos(1,k),pos(2,k),pos(3,k),'k.','MarkerSize',15)
        hold on
    end
end
xlabel('x')
ylabel('y')
view(-27,57)


