clear all
close all
%% trayectorias
tf = 20;
t0=0;
dt = 0.1;
time=t0:dt:tf;    

pos = zeros(6,length(time));
vel = zeros(6,length(time));
ace = zeros(6,length(time));

%% helicoide 3D (done) (revisada para analisis, sin interferencias)
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

%% cicloide 2D (done) (revisada para analisis, libre de interferencias y singularidades)

a = 0.25;
b = 0.9;
ini = [3 2];
inclin = pi/7;
Q = [cos(inclin) -sin(inclin); sin(inclin) cos(inclin)];

pos(1,:) = a*time(:)-b*sin(time(:));
pos(2,:) = a-b*cos(time(:));

pos(1:2,:)=Q*pos(1:2,:);
pos(1,:) = pos(1,:)+ini(1);
pos(2,:) = pos(2,:)+ini(2);

[pos(3,:),vel(3,:),ace(3,:)]=gen_coeff(0,0.3,0,0,-0.08,0,t0,tf,dt);

vel(1,:) = a - b*cos(time(:));
vel(2,:) = b*sin(time(:));
vel(1:2,:)=Q*vel(1:2,:);

ace(1,:) = b*sin(time(:));
ace(2,:) = b*cos(time(:));
ace(1:2,:)=Q*ace(1:2,:);



%% lemniscata de Bernouilli 2D (done) (revisado para analisis)

d = 2;
center = [5 4];
inclin = -pi/2;
Q = [cos(inclin) -sin(inclin); sin(inclin) cos(inclin)];

pos(1,:) = d*sqrt(2).*cos(time(1,:))./((sin(time(1,:)).^2)+1);
pos(2,:) = d*sqrt(2).*cos(time(1,:)).*sin(time(1,:))./((sin(time(1,:)).^2)+1);
pos(1:2,:)=Q*pos(1:2,:);
pos(1,:) = center(1)+pos(1,:);
pos(2,:) = center(2)+pos(2,:);

vel(1,:) = - (2^(1/2)*d.*sin(time(1,:)))./(sin(time(1,:)).^2 + 1) - (2*2^(1/2)*d.*cos(time(1,:)).^2.*sin(time(1,:)))./(sin(time(1,:)).^2 + 1).^2;
vel(2,:) = (2^(1/2)*d.*cos(time(1,:)).^2)./(sin(time(1,:)).^2 + 1) - (2^(1/2)*d.*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1) - (2*2^(1/2)*d.*cos(time(1,:)).^2.*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1).^2;
vel(1:2,:)=Q*vel(1:2,:);

ace(1,:) = (6*2^(1/2)*d.*cos(time(1,:)).*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1).^2 - (2*2^(1/2)*d.*cos(time(1,:)).^3)./(sin(time(1,:)).^2 + 1).^2 - (2^(1/2)*d.*cos(time(1,:)))./(sin(time(1,:)).^2 + 1) + (8*2^(1/2)*d.*cos(time(1,:)).^3.*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1).^3;
ace(2,:) = (6*2^(1/2)*d.*cos(time(1,:)).*sin(time(1,:)).^3)./(sin(time(1,:)).^2 + 1).^2 - (6*2^(1/2)*d.*cos(time(1,:)).^3.*sin(time(1,:)))./(sin(time(1,:)).^2 + 1).^2 + (8*2^(1/2)*d.*cos(time(1,:)).^3.*sin(time(1,:)).^3)./(sin(time(1,:)).^2 + 1).^3 - (4*2^(1/2)*d.*cos(time(1,:)).*sin(time(1,:)))./(sin(time(1,:)).^2 + 1);
ace(1:2,:)=Q*ace(1:2,:);

pos(3,:) = 0.1.*cos(time(1,:));
vel(3,:) = -0.1.*sin(time(1,:));
ace(3,:) = -0.1.*cos(time(1,:));


%% lemniscata de Bernouilli 3D (revisada para analisis)

d = 2;
center = [5 5];
inclin = -pi/6;
Q = [cos(inclin) -sin(inclin); sin(inclin) cos(inclin)];

pos(1,:) = d*sqrt(2).*cos(time(1,:))./((sin(time(1,:)).^2)+1);
pos(2,:) = d*sqrt(2).*cos(time(1,:)).*sin(time(1,:))./((sin(time(1,:)).^2)+1);
pos(1:2,:)=Q*pos(1:2,:);
pos(1,:) = center(1)+pos(1,:);
pos(2,:) = center(2)+pos(2,:);

vel(1,:) = - (2^(1/2)*d.*sin(time(1,:)))./(sin(time(1,:)).^2 + 1) - (2*2^(1/2)*d.*cos(time(1,:)).^2.*sin(time(1,:)))./(sin(time(1,:)).^2 + 1).^2;
vel(2,:) = (2^(1/2)*d.*cos(time(1,:)).^2)./(sin(time(1,:)).^2 + 1) - (2^(1/2)*d.*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1) - (2*2^(1/2)*d.*cos(time(1,:)).^2.*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1).^2;
vel(1:2,:)=Q*vel(1:2,:);

ace(1,:) = (6*2^(1/2)*d.*cos(time(1,:)).*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1).^2 - (2*2^(1/2)*d.*cos(time(1,:)).^3)./(sin(time(1,:)).^2 + 1).^2 - (2^(1/2)*d.*cos(time(1,:)))./(sin(time(1,:)).^2 + 1) + (8*2^(1/2)*d.*cos(time(1,:)).^3.*sin(time(1,:)).^2)./(sin(time(1,:)).^2 + 1).^3;
ace(2,:) = (6*2^(1/2)*d.*cos(time(1,:)).*sin(time(1,:)).^3)./(sin(time(1,:)).^2 + 1).^2 - (6*2^(1/2)*d.*cos(time(1,:)).^3.*sin(time(1,:)))./(sin(time(1,:)).^2 + 1).^2 + (8*2^(1/2)*d.*cos(time(1,:)).^3.*sin(time(1,:)).^3)./(sin(time(1,:)).^2 + 1).^3 - (4*2^(1/2)*d.*cos(time(1,:)).*sin(time(1,:)))./(sin(time(1,:)).^2 + 1);
ace(1:2,:)=Q*ace(1:2,:);

[pos(3,:),vel(3,:),ace(3,:)]=gen_coeff(1,2,0,0,0.7,0,t0,tf,dt);

[pos(4:6,:),vel(4:6,:),ace(4:6,:)] = quinticpolytraj([0 0 0; 0.1 -0.2 0; 0 -0.05 0.3; 0.3 0.1 0]',[t0 5 12 tf],time);

%% trayectoria tramos 2D (estrella) (revisada para analisis, sin colisiones)
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
times = [t0 4 8 12 16 tf];
[pos(1:2,:), vel(1:2,:), ace(1:2,:)] = quinticpolytraj(waypoint,times,time);

[pos(3,:), vel(3,:), ace(3,:)] = quinticpolytraj([0 -0.3 0.1 -0.15 0.1 0],times,time);

%% trayectoria tramos 3D (revisada para analisis, sin  colisiones)

waypoint = [5 5 6; 5 5 4; 7 7 4; 7 3 4;7 3 2; 3 3 3; 3 6 3; 5 7 3; 5 7 6; 5 5 6]';
times = [t0 2 4 6.5 8.5 11.3 13.5 16 18 tf];
[pos(1:3,:), vel(1:3,:), ace(1:3,:)] = quinticpolytraj(waypoint,times,time);
[pos(4:6,:), vel(4:6,:), ace(4:6,:)] = quinticpolytraj([0 0 0; -0.1 0 0; 0.1 0 0; 0 0.2 0; 0 0 0.05; 0 0 0.3; 0 -0.1 0; 0 0.1 0; 0.4 0 0; 0.1 0.2 0]',times,time);

%% plot 3d
figure(1)
plot3(pos(1,:),pos(2,:),pos(3,:),'k.--','MarkerEdgeColor','red')
axis([0 10 0 10 0 10])
grid on
xlabel('x')
ylabel('y')
zlabel('z')

%% plot 2d
figure(1)
plot(pos(1,:),pos(2,:),'k.--','MarkerEdgeColor','red')
axis([0 10 0 10])
grid on
xlabel('x')
ylabel('y')



