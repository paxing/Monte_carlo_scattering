%2.2.7 Photon cut-off et Roulette
function [W] = Roulette(W,m)
P=rand(1);
if P<=(1/m)
    W=W*m; 
else 
    W=0;
end

