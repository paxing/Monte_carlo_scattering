%2.2.4 Propagation
function [x,y,z] = Propagation(L,x,y,z,mu_x,mu_y,mu_z)
x=x+ceil(mu_x*L);
y=y+ceil(mu_y*L);
z=z+ceil(mu_z*L);