%2.2.5 Absorption
function [W,A] = Absorption(A,W,mu_t,mu_a)
Delta_W=(mu_a/mu_t)*W;
W=W-Delta_W;
A=A+Delta_W;
