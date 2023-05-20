%Diffusion
function [mu_x,mu_y,mu_z] = Diffusion(phi, cos_theta,mu_x,mu_y,mu_z)

if mu_z==1;
    
    mu_z=cos_theta;
    mu_x=sin(acos(cos_theta))*sin(phi);
    mu_y=sin(acos(cos_theta))*cos(phi);
else
    mu_z=-sqrt(1-mu_z^2)*sin(acos(cos_theta))*cos(phi)+mu_z*cos_theta;
    mu_x=(sin(acos(cos_theta))*(mu_x*mu_z*cos(phi)-mu_y*sin(phi)))/(sqrt(1-mu_z^2))+mu_x*cos_theta;
    mu_y=(sin(acos(cos_theta))*(mu_y*mu_z*cos(phi)-mu_x*sin(phi)))/(sqrt(1-mu_z^2))+mu_y*cos_theta;
end
