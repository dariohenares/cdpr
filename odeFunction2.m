function dz = odeFunction2(t,z,ladoderecho)

u = z(1:3); %valores de posicion
v = z(4:6); %valores de velocidad

dz = [v; ladoderecho];