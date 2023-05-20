%2.2.6 Diffusion par les tissus
%Angles
function [cos_theta,phi] = Angles(g)
P2=rand(1);
P3=rand(1);
phi=2*pi*P2;

cos_theta=(1/(2*g))*(1+g^2-((1-g^2)/(1-g+2*g*P3))^2);
