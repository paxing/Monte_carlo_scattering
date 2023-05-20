clear all, clc;


load('BlackInk.mat');
load('Proprietes.mat');
mua_data=readtable('mu_a um-1.csv');
lambda=mua_data{:,1};
mua_sc=mua_data{:,8};
mua_epidermis=mua_data{:,9};
mua_dermis=mua_data{:,12};
res_factor=8; %meaning 1 node = 1 micrometer
unit_convert=res_factor; %NOT FINAL TBD
%% Beer-Lambert en 3d 3.3




%paramètres

idx=find(lambda==694);
%Selection longueur d'onde
lbd=lambda(idx);
%************************************************************************
%absroption
%************************************************************************
mu_a_sc=mua_sc(idx)/unit_convert;
mu_a_ed=mua_epidermis(idx)/unit_convert;
mu_a_d=mua_dermis(idx)/unit_convert;
mu_a_tatou=Blackink.absorption(idx)/unit_convert;

%************************************************************************
%scattering en noeuds
%************************************************************************
mu_s_sc=52/unit_convert;
mu_s_ed=90/unit_convert;
mu_s_d=60/unit_convert;
mu_s_tatou=90/unit_convert;

%%
% initialisation paramètre de simulation
N=100000;

W=1;
W_c=0.01;
m=10;

tatou=true;
if tatou==false
    % Dimensions des couches en noeuds selon z
    N_sc=12*res_factor;
    N_ed=80*res_factor;
    N_d=500*res_factor;
    
    %taille de maillage en z
    Nz=N_sc+N_ed+N_d;
    Nx=600;
    Ny=600;
    A=zeros(Nx,Ny,Nz);
    
    
    
    %abs en noeuds
    % mu_a_sc=4.733/um_to_nodes;
    % mu_a_ed=40/um_to_nodes;
    % mu_a_d=4/um_to_nodes;
    
    
    
    %indice de réfraction
    
    n_sc=1.42;
    n_ed=1.42;
    n_d=1.39;
    
    %identifcation des propriétés selon le maillage
    mu_at=[mu_a_sc*ones(N_sc,1);mu_a_ed*ones(N_ed,1);mu_a_d*ones(N_d,1)];
    mu_st=[mu_s_sc*ones(N_sc,1);mu_s_ed*ones(N_ed,1);mu_s_d*ones(N_d,1)];
    n_t=[n_sc*ones(N_sc,1);n_ed*ones(N_ed,1);n_d*ones(N_d,1)];
    g=0.7;
elseif tatou==true
    N_sc=12*res_factor;
    N_ed=80*res_factor;
    N_d1=235*res_factor;
    N_tatou=30*res_factor;
    N_d2=235*res_factor;
    
    %taille de maillage en z
    Nz=N_sc+N_ed+N_d1+N_tatou+N_d2;
    Nx=600;
    Ny=600;
    A=zeros(Nx,Ny,Nz);
    
    mu_at=[mu_a_sc*ones(N_sc,1);mu_a_ed*ones(N_ed,1);mu_a_d*ones(N_d1,1);mu_a_tatou*ones(N_tatou,1);mu_a_d*ones(N_d2,1)];
    mu_st=[mu_s_sc*ones(N_sc,1);mu_s_ed*ones(N_ed,1);mu_s_d*ones(N_d1,1);mu_s_tatou*ones(N_tatou,1);mu_s_d*ones(N_d2,1)];
    g=0.7;
    
end
%%


[A]=Monte_carlo(A,mu_at,mu_st,g,Nx,Ny,Nz,N,W_c,m);


%% test



%%
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


%%
figplot=figure()
conc=sum(A,2);
conc=reshape(conc,Ny,Nz);
imagesc(linspace(0,Ny,Ny),linspace(0,Nz,Nz),log10(transpose(conc)))
cb=colorbar
cb.Label.String = 'log(Fr\''{e}quence)';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
hold on
line([0, Nx], [N_sc,N_sc], 'Color', 'w');
hold on
line([0, Nx], [N_ed+N_sc,N_ed+N_sc], 'Color', 'w');
hold on
if tatou==true
    line([0, Nx], [N_ed+N_sc+N_d1,N_ed+N_sc+N_d1], 'Color', 'w');
    hold on
    line([0, Nx], [N_ed+N_sc+N_d1+N_tatou,N_ed+N_sc+N_d1+N_tatou], 'Color', 'w');
end
xlabel('Position en x $$  (\mu m) $$','Interpreter','latex')
ylabel('Position en z $$  (\mu m) $$','Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',28)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
%print(figplot,'test_3couche','-dpng','-r0')

%%
figplot=figure()
Az=sum(conc,1);

plot(linspace(0,2,M),Az,'LineWidth',3)
xlabel('Position en z $$  (mm) $$','Interpreter','latex')
ylabel('Fr\''{e}quence ','Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',28)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
%print(figplot,'diff_3d_z','-dpng','-r0')
%% Modèle à 7 couches

clear all, clc;
% ALL UNITS MUST BE IN MICROMETER (default 1 NODE = 1 MICROMETER)

res_factor=2; %meaning 1 node = 0.5 micrometer
unit_convert=10000*res_factor; %NOT FINAL

%********************************
%EPIDERMIS
%********************************
%Stratum corneum
sc_d=15*res_factor;
sc_mua=4/unit_convert;
sc_mus=10/unit_convert;
sc_g=0.8;

%Stratum granulosum
sl_d=5*res_factor;
sl_mua=50/unit_convert;
sl_mus=2/unit_convert;
sl_g=0.75;

%Stratum malpighi 
sm_d=150*res_factor;
sm_mua=2/unit_convert;
sm_mus=2/unit_convert;
sm_g=0.7;

%stratum basale
sb_d=10*res_factor;
sb_mua=1/unit_convert;
sb_mus=9/unit_convert;
sb_g=0.87;

%********************************
%DERMIS
%********************************
dermis_d=500*res_factor;
dermis_mua=0.4/unit_convert;
dermis_mus=15/unit_convert;
dermis_g=0.9;

%********************************
%HYPODERMIS
%********************************

va_d=200*res_factor;
va_mua=50/unit_convert;
va_mus=3/unit_convert;
va_g=0.9;

%adipous
ad_d=400*res_factor;
ad_mua=1/unit_convert;
ad_mus=15/unit_convert;
ad_g=0.6;


%********************************
%TATTOO
%********************************

tat_d=30*res_factor;
tat_mua=100/unit_convert;
tat_mus=10/unit_convert;
tat_g=0.6;


%% skin initialization

%initalization of absorption matrix
Nz=sc_d+sl_d+sm_d+sb_d+dermis_d+va_d+ad_d;
Nx=300;
Ny=300;
A=zeros(Nx,Ny,Nz);

%proprties in matrix form
skin=[sc_d, sl_d, sm_d, sb_d, dermis_d, va_d, ad_d;
    sc_mua, sl_mua, sm_mua, sb_mua, dermis_mua, va_mua, ad_mua;
    sc_mus, sl_mus, sm_mus, sb_mus, dermis_mus, va_mus, ad_mus;
    sc_g, sl_g, sm_g, sb_g, dermis_g, va_g, ad_g];


mu_at=[];
mu_st=[];
gt=[];

for i=1:7
    mu_at=[mu_at;skin(2,i)*ones(skin(1,i),1)];
    mu_st=[mu_st;skin(3,i)*ones(skin(1,i),1)];
    gt=[gt;skin(4,i)*ones(skin(1,i),1)];
end

%switch for the tattoo
tattoo=false;
if tattoo
    %calculation of the position of the tattoo in nodes
    tat_up=sc_d+sl_d+sm_d+sb_d+round(dermis_d/2)-round(tat_d/2);
    tat_down=sc_d+sl_d+sm_d+sb_d+round(dermis_d/2)+round(tat_d/2);
    
    %change parameter values 
    mu_at(tat_up:tat_down)=tat_mua;
    mu_st(tat_up:tat_down)=tat_mus;
    gt(tat_up:tat_down)=tat_g;
end

%% simulation

% initialisation paramètre de simulation
N=1000000;
W_c=0.01;
m=10;

[A]=Monte_carlo(A,mu_at,mu_st,gt,Nx,Ny,Nz,N,W_c,m);




%% plot 2d figure
figplot=figure()
conc=mean(A,2);
%conc=A(round(Nx/2),:,:);
conc=reshape(conc,Ny,Nz);
imagesc(linspace(0,Ny,Ny)/res_factor,linspace(0,Nz,Nz)/res_factor,log10(transpose(abs(conc))))
cb=colorbar
colormap jet
cb.Label.String = 'log(Fr\''{e}quence)';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

xlabel('Position en x $$  (\mu m) $$','Interpreter','latex')
ylabel('Position en z $$  (\mu m) $$','Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',28)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
%print(figplot,'test_3couche','-dpng','-r0')

%%

%% plot 3d figures



%conver to log10 value 
intensity=log10(abs(A(round(Nx/2):end,:,:)));



%return the index of each value in A

idx=find(A(round(Nx/2):end,:,:));
%give the x,y,z position for each node
[X, Y, Z] = ind2sub(size(A(round(Nx/2):end,:,:)), idx);



f=figure()
pointsize = 2;

scatter3(X(:), Y(:), Z(:), pointsize, intensity(idx));
cb=colorbar
colormap jet
cb.Label.String = 'log(Fr\''{e}quence)';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

set(gca, 'ZDir','reverse')
ax=gca;
ax.LineWidth=2;
xlabel('Position en x $$  (\mu m) $$','Interpreter','latex')
set(get(gca,'xlabel'),'rotation',10)

zlabel('Position en z $$  (\mu m) $$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',18)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
%caxis([-6 0.5])
%print(f,'test_3couches_tat','-dpng','-r0')

%%
