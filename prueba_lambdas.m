close all

%% generacion de trayectorias
tf = 20;
t0=0;
dt = 0.1;
time=t0:dt:tf;    

pos = zeros(6,length(time));
vel = zeros(6,length(time));
ace = zeros(6,length(time));

r = 2.5;
omega = 1;
center = [5 5 1];
v_0 = 0.3;

pos(1,:) = center(1)+r*cos(omega*time(:));
pos(2,:) = center(2)+r*sin(omega*time(:));
pos(3,:) = center(3)+v_0*time(:);

pos(4:6,:) = zeros(3,length(time));
vel(4:6,:) = zeros(3,length(time));
ace(4:6,:) = zeros(3,length(time));

vel(1,:) = -omega*r*sin(omega*time(:));
vel(2,:) = r*omega*cos(omega*time(:));
vel(3,:) = v_0;

ace(1,:) = -omega^2*r*cos(omega*time(:));
ace(2,:) = -r*omega^2*sin(omega*time(:));
ace(3,:) = 0;

[pos(4:6,:), vel(4:6,:), ace(4:6,:)] = quinticpolytraj([0 0 0; 0.2 -0.15 0.15; 0 0 0]',[t0 10 tf],time);

%% definicion de parametros invariables e inicializacion

V = [-0.5 -0.5 -0.5; -0.5 0.5 -0.5; 0.5 0.5 -0.5 ; 0.5 -0.5 -0.5;...
    -0.5 -0.5 0.5; -0.5 0.5 0.5; 0.5 0.5 0.5 ; 0.5 -0.5 0.5]';
A = [0 0 0; 0 10 0; 10 10 0; 10 0 0; 0 0 10; 0 10 10; 10 10 10; 10 0 10]';
B = [-0.5 -0.05 0.5; -0.5 0.05 0.5; 0.5 0.05 0.5; 0.5 -0.05 0.5;...
    -0.05 -0.5 -0.5; -0.05 0.5 -0.5; 0.05 0.5 -0.5; 0.05 -0.5 -0.5]';

t_min = 0;
t_max = 1000;
lb = ones(8,1)*t_min;
ub = ones(8,1)*t_max;


m = 100; 
I_p = [16.67 0 0; 0 16.67 0; 0 0 16.67]; 
f_g = [0; 0; -m*9.81; 0; 0; 0];

t_linprog = zeros(8,length(time));
t_quadprog = zeros(8,length(time));
t_penrose = zeros(8,length(time));

pos_ode = zeros(6,length(time)); 
vel_ode = zeros(6,length(time)); 
ace_ode = zeros(6,length(time));

pos_ode(:,1) = pos(:,1);
vel_ode(:,1) = vel(:,1);
ace_ode(:,1) = ace(:,1);

dist = zeros(4,length(time));
dist_p = zeros(16,length(time));

cubo_plot = [];
cables_plot = [];
%% Blucle de cálculo - progagación en tiempo

for i=1:length(time)

    Q = [cos(pos(5,i))*cos(pos(6,i)) -cos(pos(5,i))*sin(pos(6,i)) sin(pos(5,i));
        cos(pos(4,i))*sin(pos(6,i))+sin(pos(4,i))*sin(pos(5,i))*cos(pos(6,i))...
        cos(pos(4,i))*cos(pos(6,i))-sin(pos(4,i))*sin(pos(5,i))*sin(pos(6,i))...
        -sin(pos(4,i))*cos(pos(5,i));
        sin(pos(4,i))*sin(pos(6,i))-cos(pos(4,i))*sin(pos(5,i))*cos(pos(6,i))...
        sin(pos(4,i))*cos(pos(6,i))+cos(pos(4,i))*sin(pos(5,i))*sin(pos(6,i))...
        cos(pos(4,i))*cos(pos(5,i))];    
    
    l = repmat(pos(1:3,i),1,8)+Q*B-A; 
    L = [norm(l(:,1)) norm(l(:,2)) norm(l(:,3)) norm(l(:,4))...
         norm(l(:,5)) norm(l(:,6)) norm(l(:,7)) norm(l(:,8))]; 
    J = [(l(1:3,:)./L)' (cross(Q*B,repmat(pos(1:3,i),1,8)-A)./L)'];
    W = -J'; 
    K = Q*I_p*Q'; 
    M = [m*eye(3) zeros(3,3); zeros(3,3) K];
    C = [zeros(3,1); cross(vel(4:6,i),K*vel(4:6,i))]; 
    
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

    t_pinv = pinv(W)*(M*ace(:,i) + C - f_g);
    t_penrose(:,i) = t_pinv;
    N = null(W);
    
    lambda_l = optimvar('lambda_l',2);
    prob_l = optimproblem('Objective',ones(1,8)*(t_pinv+N*lambda_l),'ObjectiveSense','min');
    prob_l.Constraints.c1 = W*(t_pinv+N*lambda_l) == M*ace(:,i)+C-f_g;
    prob_l.Constraints.c2 = ones(8,1)*t_min-t_pinv <= N*lambda_l;
    prob_l.Constraints.c3 = ones(8,1)*t_max-t_pinv >= N*lambda_l;
    problem_l = prob2struct(prob_l);
    lambda_l = linprog(problem_l);
    if isempty(lambda_l) == 1
        t_linprog(:,i) = zeros(8,1);
    else
        t_linprog(:,i) = t_pinv+N*lambda_l;
    end

%     lambda_q = optimvar('lambda_q',2);
%     prob_q = optimproblem('Objective',(N*lambda_q)'*(N*lambda_q),'ObjectiveSense','min');
%     prob_q.Constraints.c1 = W*(t_pinv+N*lambda_q) == M*ace(:,i)+C-f_g;
%     prob_q.Constraints.c2 = ones(8,1)*t_min-t_pinv <= N*lambda_q;
%     prob_q.Constraints.c3 = ones(8,1)*t_max-t_pinv >= N*lambda_q;
%     problem_q = prob2struct(prob_q);
%     problem_q.options.LinearSolver = 'dense';
%     problem_q.options.MaxIterations = 1000;
%     problem_q.options.ConstraintTolerance = 1e-4;  
%     lambda_q = quadprog(problem_q);
%     if isempty(lambda_q) == 1
%         t_quadprog(:,i) = zeros(8,1);
%     else
%         t_quadprog(:,i) = t_pinv+N*lambda_q;
%     end
    
    %programacion cuadratica
    
    Aeq = W; 
    Beq = M*ace(:,i)+C-f_g;
    t_q = quadprog(eye(8),zeros(8,1),[],[],Aeq,Beq,lb,ub); %t_opt  vector de tensiones de los cables
    if isempty(t_q) == 1
        t_quadprog(:,i) = zeros(8,1);
    else
        t_quadprog(:,i)= t_q; 
    end




%     %parametros para plotear en tiempo real
% 
%     cubo = [(pos(1:3,i)+Q*V(:,1))' (pos(1:3,i)+Q*V(:,2))';
%             (pos(1:3,i)+Q*V(:,2))' (pos(1:3,i)+Q*V(:,3))';
%             (pos(1:3,i)+Q*V(:,3))' (pos(1:3,i)+Q*V(:,4))';
%             (pos(1:3,i)+Q*V(:,4))' (pos(1:3,i)+Q*V(:,1))';
%             (pos(1:3,i)+Q*V(:,5))' (pos(1:3,i)+Q*V(:,6))';
%             (pos(1:3,i)+Q*V(:,6))' (pos(1:3,i)+Q*V(:,7))';
%             (pos(1:3,i)+Q*V(:,7))' (pos(1:3,i)+Q*V(:,8))';
%             (pos(1:3,i)+Q*V(:,8))' (pos(1:3,i)+Q*V(:,5))';
%             (pos(1:3,i)+Q*V(:,1))' (pos(1:3,i)+Q*V(:,5))';
%             (pos(1:3,i)+Q*V(:,2))' (pos(1:3,i)+Q*V(:,6))';
%             (pos(1:3,i)+Q*V(:,3))' (pos(1:3,i)+Q*V(:,7))';
%             (pos(1:3,i)+Q*V(:,4))' (pos(1:3,i)+Q*V(:,8))'];
%           
%     cubo_plot = cat(2,cubo_plot,cubo);
%         
%     cables = [(pos(1:3,i)+Q*B(:,1:8)); A(:,1:8)]';
%     
%     cables_plot = cat(2,cables_plot,cables);
        
%     %% ODE45
%     % cambiar t_linprog, t_quadprog, t_pinv segun convenga
%     ladoderecho = M\(W*t_l+f_g-C); 
%     
%     if i<length(time)    %para que no calcule en el instante final
%         [t_sol,z] = ode45(@(t,z) odeFunction(t,z,ladoderecho), [time(i) time(i+1)], [pos_ode(:,i); vel_ode(:,i)]);
% 
%         pos_ode(:,i+1) = z(length(z),1:6); %cojo el ultimo valor, que es el que equivale al instante i+1
%         vel_ode(:,i+1) = z(length(z),7:12);
%         ace_ode(:,i+1) = ladoderecho;
%     end

end

% %% PLOTS
% cont = 0;
% figure(1)
% set(gcf, 'Position',  [100, 100, 1400, 600])
% for i=1:6:6*length(time)
%     subplot(1,2,1)
%     for k = 1:12
%         plot3(cubo_plot(k,i:3:i+3),cubo_plot(k,i+1:3:i+4),cubo_plot(k,i+2:3:i+5),'k');
%         hold on
%     end
% 
%     for j = 1:8
%         plot3(cables_plot(j,i:3:i+3),cables_plot(j,i+1:3:i+4),cables_plot(j,i+2:3:i+5))
%         hold on
%     end
%    
%     axis([0 10 0 10 0 10])
%     hold off
%     grid off
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     view (-65,18)
%     
%     %base fija
%     line([0 0],[0 10],[0 0],'color','k')
%     line([0 0],[0 0],[0 10],'color','k')
%     line([0 0],[10 10],[0 10],'color','k')
%     line([0 0],[0 10],[10 10],'color','k')
%     line([10 10],[0 10],[0 0],'color','k')
%     line([10 10],[0 0],[0 10],'color','k')
%     line([10 10],[10 10],[0 10],'color','k')
%     line([10 10],[0 10],[10 10],'color','k')
%     line([0 10],[10 10],[10 10],'color','k')
%     line([0 10],[0 0],[10 10],'color','k')
%     line([0 10],[0 0],[0 0],'color','k')
%     line([0 10],[10 10],[0 0],'color','k')
%     
%     cont=cont+1;
%     
%     txt1=['Distancia minima entre cable 1 y 5: ' num2str(dist(1,cont)) ' m '];
%     txt2=['Distancia minima entre cable 2 y 6: ' num2str(dist(2,cont)) ' m '];
%     txt3=['Distancia minima entre cable 3 y 7: ' num2str(dist(3,cont)) ' m '];
%     txt4=['Distancia minima entre cable 4 y 8: ' num2str(dist(4,cont)) ' m '];
%     text1_p=['cable 1 - arista opuesta: ',num2str(dist_p(1,cont)),' m ; cable 1 - arista lateral: ',num2str(dist_p(2,cont)),' m'];
%     text2_p=['cable 2 - arista opuesta: ',num2str(dist_p(3,cont)),' m ; cable 2 - arista lateral: ',num2str(dist_p(4,cont)),' m'];
%     text3_p=['cable 3 - arista opuesta: ',num2str(dist_p(5,cont)),' m ; cable 3 - arista lateral: ',num2str(dist_p(6,cont)),' m'];
%     text4_p=['cable 4 - arista opuesta: ',num2str(dist_p(7,cont)),' m ; cable 4 - arista lateral: ',num2str(dist_p(8,cont)),' m'];
%     text5_p=['cable 5 - arista opuesta: ',num2str(dist_p(9,cont)),' m ; cable 5 - arista lateral: ',num2str(dist_p(10,cont)),' m'];
%     text6_p=['cable 6 - arista opuesta: ',num2str(dist_p(11,cont)),' m ; cable 6 - arista lateral: ',num2str(dist_p(12,cont)),' m'];
%     text7_p=['cable 7 - arista opuesta: ',num2str(dist_p(13,cont)),' m ; cable 7 - arista lateral: ',num2str(dist_p(14,cont)),' m'];
%     text8_p=['cable 8 - arista opuesta: ',num2str(dist_p(15,cont)),' m ; cable 8 - arista lateral: ',num2str(dist_p(16,cont)),' m'];
%     hold off
%     
%     subplot(1,2,2)
%     plot(time(1:cont),t_quadprog(:,1:cont))
%     text(1,1600,{txt1 txt2 txt3 txt4 text1_p text2_p text3_p text4_p text5_p text6_p text7_p text8_p},'FontSize',8)
%     axis([0 20 0 2000])
%     hold off
%     grid on
%     xlabel('Tiempo (s)')
%     ylabel('Tensiones (N)')
%     legend('1','2','3','4','5','6','7','8')
%     pause(0.05)
% end


