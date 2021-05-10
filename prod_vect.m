function [vect_z]=prod_vect(v1,v2)
%para vectores en el plano (Z=0)
%v1,v2 dim(2,1)
vect_z = v1(1)*v2(2)-v1(2)*v2(1);


