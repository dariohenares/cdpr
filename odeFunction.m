function dz = odeFunction(t,z,ladoderecho)

u = z(1:6); %valores de posicion
v = z(7:12); %valores de velocidad

dz = [v; ladoderecho]; %lado derecho se calcula en el principal para evitar tener que llamar a muchas variables a la funcion



