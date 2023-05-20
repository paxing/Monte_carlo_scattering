%% Modèle à 7 couches

clear all, clc;

% ALL UNITS MUST BE IN MICROMETER (default 1 NODE = 1 MICROMETER)
load('BlackInk.mat');
load('Proprietes.mat');
mua_data=readtable('mu_a um-1_old.csv');
lambda=mua_data{:,1};
mua_sc=mua_data{:,8};
mua_epidermis=mua_data{:,9};
mua_dermis=mua_data{:,12};
res_factor=8; %meaning 1 node = 1 micrometer
unit_convert=res_factor; %NOT FINAL TBD

%****************************************************
%EPIDERMIS
%****************************************************
%Stratum corneum
sc_d=15*res_factor;
sc_mua=4/unit_convert;
sc_mus=10/unit_convert;
sc_g=0.8;

%Stratum granulosum
sl_d=5*res_factor;
sl_mua=4/unit_convert;
sl_mus=5/unit_convert;
sl_g=0.75;

%Stratum malpighi 
sm_d=150*res_factor;
sm_mua=2/unit_convert;
sm_mus=12/unit_convert;
sm_g=0.87;

%stratum basale
sb_d=10*res_factor;
sb_mua=1/unit_convert;
sb_mus=9/unit_convert;
sb_g=0.92;

%****************************************************
%DERMIS
%****************************************************
dermis_d=300*res_factor;
dermis_mua=2/unit_convert;
dermis_mus=15/unit_convert;
dermis_g=0.9;

%****************************************************
%HYPODERMIS
%****************************************************

va_d=100*res_factor;
va_mua=10/unit_convert;
va_mus=13/unit_convert;
va_g=0.89;

%adipeous
ad_d=500*res_factor;
ad_mua=2/unit_convert;
ad_mus=15/unit_convert;
ad_g=0.86;


%****************************************************
%TATTOO
%****************************************************
tat_d=30*res_factor;
tat_mua=120/unit_convert;
tat_mus=10/unit_convert;
tat_g=0.85;


%% skin initialization

%initalization of absorption matrix
Nz=sc_d+sl_d+sm_d+sb_d+dermis_d+va_d+ad_d;
Nx=250;
Ny=250;
A=zeros(Nx,Ny,Nz); %initialize all absorption event at 0

%*************************************************
%propreties of skin in matrix form (only 1D in z axis)
skin=[sc_d, sl_d, sm_d, sb_d, dermis_d, va_d, ad_d;
    sc_mua, sl_mua, sm_mua, sb_mua, dermis_mua, va_mua, ad_mua;
    sc_mus, sl_mus, sm_mus, sb_mus, dermis_mus, va_mus, ad_mus;
    sc_g, sl_g, sm_g, sb_g, dermis_g, va_g, ad_g];
%*************************************************



%**********************************************
%7 layers model in 1D (z axis)
%**********************************************
mu_at=[];
mu_st=[];
gt=[];

for i=1:7
    mu_at=[mu_at;skin(2,i)*ones(skin(1,i),1)];
    mu_st=[mu_st;skin(3,i)*ones(skin(1,i),1)];
    gt=[gt;skin(4,i)*ones(skin(1,i),1)];
end

%switch for the tattoo
tattoo=true;

%this add a tattou in the middle of the dermis
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

%note: new version of the Monte_carlo (add variation in xy plane)
[A]=Monte_carlo2(A,mu_at,mu_st,gt,Nx,Ny,Nz,N,W_c,m,5);




%% plot 2d figure

figplot=figure()
%this calculate the mean value in the xz plane
conc=mean(A,2);
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
print(f,'test_7couches_tat','-dpng','-r0')