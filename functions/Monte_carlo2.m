function [A]=Monte_carlo2(A,mu_at,mu_st,gt,Nx,Ny,Nz,N,W_c,m,period)


for i=1:N
    
    %initialisation du photon
    z=0;
    y=ceil(Ny/2);
    x=ceil(Nx/2);
    W=1;
    mu_z=1;
    mu_y=0;
    mu_x=0;
    mu_a=mu_at(1);
    mu_s=mu_st(1);
    mu_t=mu_a+mu_s;
    g=gt(1);

    while W ~=0
        
        %2.2.3 Distance de parcours
        L=((-log(rand(1))/mu_t)); %distance en noeuds
        %2.2.4 propagation du photon
        
        [x,y,z] = Propagation(L,x,y,z,mu_x,mu_y,mu_z);
       
               
        %2.2.5 absorption du photon
        if z<=Nz && y<=Ny && x<=Nx && z>0 && y>0 && x> 0
            mu_a=mu_at(z);
            mu_s=mu_st(z);
            mu_t=mu_a+mu_s;
            g=gt(z);
            
            %add of variation in xy plane
            mu_a=mu_a*(1+1/10*sin(period/Nx*x)*sin(period/Ny*y));
            mu_s=mu_s*(1+1/10*sin(period/Nx*x)*sin(period/Ny*y));
            g=g*(1+1/20*sin(period/Nx*x)*sin(period/Ny*y));
            if g>1
                g=1;
            end
            
            [W,A(x,y,z)] = Absorption(A(x,y,z),W,mu_t,mu_a);
        else
            W=0;
        end 
        
        
        %2.2.6 dffusion
        [cos_theta,phi] = Angles(g);%calcul angle
        
        [mu_x,mu_y,mu_z] = Diffusion(phi,cos_theta,mu_x,mu_y,mu_z);
        
        %2.2.7 Photon cut-off et Roulette
        if W < W_c % on part la roulette
            W=Roulette(W,m);
       
        end
    end
    
    
end