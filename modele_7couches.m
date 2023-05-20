%%
close all;clear; clc;
lambda = (400:1:1100)'; %nm
g = 0.8*ones(length(lambda), 1);
%% Stratum corneum
% Faire une courbe à partir des données de Tuchin et Meglinski & Matcher
% (2000).
% On utilise seulement les deux points connus dans le visible. On n'utilise
% pas les points dans l'UV comme on sait que les propriétés d'absorption et
% de scattering sont fort différentes. À partir de ces deux points, on fit
% une rationnelle comme on sait que les courbes obéissent ainsi, d'après
% les data de l'épithélium et du derme.

lambda_2 = [400 633];
mus_sc2 = [2000 1000]*0.0001;

% fit obtenu : f(x) = 46.44 / (x -168.1)

mus_sc = 46.44 ./ (lambda - 168.1);

figure(1)
scatter(lambda_2, mus_sc2, 'filled')
ylabel('µ_s [µm^{-1}]')
xlabel('Longueur d''onde [nm]')
hold on
plot(lambda, mus_sc, 'Linewidth', 2)
legend('Data', 'Fit')

Proprietes.stratum.scattering = mus_sc;
Proprietes.stratum.g = g;
Proprietes.stratum.lambda = lambda;
%% Epiderme
% Hypothèse : Mu_s est constant sur toutes les sous-couches de l'épiderme. 
% Source : Tuchin (2000)
lambda_e_1 = [415 488 514 585 633 800];
mus_e_1 = [800 600 600 470 450 420]*0.0001; %um^-1

lambda_e_2 = [370 420 470 488 514 520 570 620 633 670 720 770 820 830 870 920 970 1020 1064 1070 1120];
mus_e_2 = [11.56 9.82 7.96 7.41 6.67 6.51 5.52 4.90 4.76 4.48 4.11 3.79 3.60 3.56 3.41 3.32 3.15 3.02 2.97 2.97 2.86]*0.0001; %um^-1

figure(2)
hold on
scatter(lambda_e_1, mus_e_1, 'filled')
legend()

%figure(3)
%scatter(lambda_e_2, mus_e_2, 'filled')

% Avec un g constant, les deux courbes adoptent les mêmes allures. On va
% ainsi faire un fit sur la courbe avec les vraies valeurs de mu_s pour
% avoir mu_s sur le spectre de 400 nm à 1100 nm.
% fit nous donne une rationnelle f(x) = 0.03127*x-5.551/(x - 322.4)

mus_epi = (0.03127.*lambda - 5.551)./(lambda - 322.4);

figure(2)
plot(lambda, mus_epi, 'Linewidth', 2)
ylabel('µ_s [µm^{-1}]')
xlabel('Longueur d''onde [nm]')
legend('Données (Tuchin)', 'Fit')

Proprietes.epiderme.scattering = mus_epi;
Proprietes.epiderme.g = g;
Proprietes.epiderme.lambda = lambda;
%% Derme
% Hypothèse : mu_s est le même sur toutes les couches du derme, étant donné
% l'ordre de grandeur de mu_s qui est le même dans l'article de Meglinski &
% Matcher (2000) et les incertitudes sur les courbes de Simpson & al
% (1998).
% Source: Tuchin (2000)
lambda_d_1 = [415 488 585 633 800];
mus_d_1 = [320 250 196 187.5 175]*0.0001;

figure(4)
hold on
scatter(lambda_d_1, mus_d_1, 'filled')
ylabel('µ_s [µm^{-1}]')
xlabel('Longueur d''onde [nm]')

lambda_d_2 = [420 470 488 514 520 570 620 633 670 720 770 820 830 870 920 970 1020 1064 1070 1120];
mus_d_2 = [6.85 5.36 4.90 4.32 4.2 3.5 3.07 2.99 2.78 2.54 2.33 2.18 2.15 2.05 1.99 1.90 1.84 1.70 1.79 1.74]*15./(1-0.8)*0.0001;

%figure(4)
%scatter(lambda_d_2, mus_d_2, 'filled')
%legend()

%même raisonnement que pour l'épiderme. Le fit est de
% f(x) = 8.092e07 *x^(-3.702) + 0.01556

mus_derme = 8.092e7*lambda.^(-3.702) + 0.01556;
figure(4)
plot(lambda, mus_derme, 'Linewidth', 2)
legend('Data (Tuchin)', 'Fit')

Proprietes.derme.scattering = mus_derme;
Proprietes.derme.g = g;
Proprietes.derme.lambda = lambda;

%% subcutaneous fat
% On n'a pas de données de mu_s pour le subcutaneous fat. On a seulement
% mu_s = 5 mm^-1 à lambda = 633~nm. Comme on suppose que g est constant et
% que la relation entre mu_s et mu_s' est linéaire, on va utiliser cette
% valeur pour adapter la courbe de données. À 633 nm, mu_s = 50 cm^-1, et
% mu_s' = 1.88 cm^-1, il faut donc multipler mu_s' par 50/1.88 pour trouver
% mu_s. Source : Tuchin (2000) et Meglinsi & Matcher (2002)

lambda_fat = [420 470 488 514 520 570 620 633 670 720 770 820 830 780 920 970 1020 1064 1070 1120];
mu_s_f = 50/1.88 * [4.21 3.38 3.13 2.80 2.74 2.35 1.95 1.88 1.71 1.52 1.35 1.24 1.22 1.16 1.09 1.02 0.94 0.88 0.88 0.85]*0.0001;

figure(5)
scatter(lambda_fat, mu_s_f, 'filled')
ylabel('µ_s [µm^{-1}]')
xlabel('Longueur d''onde [nm]')

%fit obtenu est de f(x) = (-0.0002774*x + 2.172)/(x - 239.1)

mus_fat = (-0.0002774.*lambda + 2.172) ./ (lambda - 239.1);

figure(5)
hold on
plot(lambda, mus_fat, 'Linewidth', 2)

Proprietes.subfat.scattering = mus_fat;
Proprietes.subfat.g = g;
Proprietes.subfat.lambda = lambda;

%%

figure(6)
hold on
plot(Proprietes.stratum.lambda, Proprietes.stratum.scattering, 'Linewidth', 2)
plot(Proprietes.epiderme.lambda, Proprietes.epiderme.scattering, 'Linewidth', 2)
plot(Proprietes.derme.lambda, Proprietes.derme.scattering, 'Linewidth', 2)
plot(Proprietes.subfat.lambda, Proprietes.subfat.scattering', 'Linewidth', 2)
ylabel('Coefficient de diffusion \mu_s [\mum^{-1}]', 'Fontsize', 12)
xlabel('Longueur d''onde [nm]', 'Fontsize', 12)
legend('Stratum corneum', 'Epidermis', 'Dermis', 'Subcutaneous Fat')

%% Tatou
% Coefficient d'absorption
load ink

lambda_black = blackink(:,1);
mu_black = blackink(:,2)*0.0001;


% curve fit black ink
% f(x) = 0.1315 /(x - 9.801)
mua_black = 0.1315 ./ (lambda - 9.801);

figure(7)
hold on
scatter(lambda_black, mu_black, 'filled')
plot(lambda, mua_black, 'Linewidth', 2)
legend('Data', 'Fit')
ylabel('Coefficient d''absorption \mua [\mum^{-1}]')
xlabel('Longueur d''onde [nm]')

Blackink.absorption = mua_black;
Blackink.scattering = 90*0.0001*ones(length(lambda), 1);
Blackink.g = 0.7*ones(length(lambda), 1);
Blackink.lambda = lambda;