clc
clearvars


% Definiera parametrar
m1 = 465; % massa för chassidelen (kg)
m2 = 55; % massa för hjulet (kg)
k1 = 5350; % fjäderkonstant för chassifjädern (N/m)
k2 = 136100; % fjäderkonstant för däckfjädern (N/m)
c1 = 310; % dämpningskonstant för chassifjädern (Ns/m)
c2 = 1250; % dämpningskonstant för däckfjädern (Ns/m)
v0 = [0; 0; 0; 0]; % initialvillkor

H = 0.27; % höjd på gupp (m)
L = 1.1; % längd på gupp (m)
v_speed = 63/3.6; % hastighet av fordonet (m/s)




% Skapa systemmatrisen A
A = [0, 1, 0, 0;
    -k1/m1, -c1/m1, 0, 0;
    0, 0, 0, 1;
    k1/m2, c1/m2, -k2/m2, -c2/m2];

% Beräkna egenvärden
eigenvalues = eig(A);

% Beräkna dt_max för varje egenvärde
dt_max_values = -2 * real(eigenvalues) ./ abs(eigenvalues).^2;

% Filtrera ut negativa värden eftersom de inte är fysiskt relevanta
dt_max_values = dt_max_values(dt_max_values > 0);

% Hitta det minsta positiva värdet
Delta_t_max = min(dt_max_values);

% Jämför med det givna tidssteget
Delta_t_given = 5e-3;
stability = Delta_t_given <= Delta_t_max;

% Visa resultat
fprintf('Det maximala stabila tidssteget Delta_t_max är: %.5f\n', Delta_t_max);
if stability
    fprintf('Det givna tidssteget %.1e uppfyller stabilitetsvillkoret.\n', Delta_t_given);
else
    fprintf('Det givna tidssteget %.1e uppfyller inte stabilitetsvillkoret.\n', Delta_t_given);
end









