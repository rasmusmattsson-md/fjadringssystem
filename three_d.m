%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created with ♡ by Hampus & Rasmus
%% Available under the MIT-license
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clearvars



% Definiera parametrar
m1 = 465; % massa för chassidelen (kg)
m2 = 55; % massa för hjulet (kg)
k1 = 5350; % fjäderkonstant för chassifjädern (N/m)
k2 = 136100 * 100; % fjäderkonstant för däckfjädern (N/m)
c1 = 310; % dämpningskonstant för chassifjädern (Ns/m)
c2 = 1250; % dämpningskonstant för däckfjädern (Ns/m)
H = 0.27; % höjd på gupp (m)
L = 1.1; % längd på gupp (m)
v_speed = 63/3.6; % hastighet av fordonet (m/s)
v0 = [0; 0; 0; 0]; % initialvillkor

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



% Initialiserar tidsvektorn och tillståndsvektorn
t_end = 1; % Definera hur lång tid simuleringen ska köra
num_steps = ceil(t_end / Delta_t_max); % Beräkna antalet tidssteg
v = zeros(4, num_steps); % Initiera tillståndsvektor för alla tidssteg
v(:, 1) = v0; % Sätt initialt tillstånd
t_span = linspace(0, t_end, num_steps); % Skapa tidsvektorn

% Euler-framåt-loop
for i = 1:num_steps-1
    % Din suspension_system funktion behöver definieras på rätt sätt för att returnera
    % derivatan av systemtillståndet. Antag att det är gjort.
    v(:, i+1) = v(:, i) + Delta_t_max * suspension_system(t_span(i), v(:, i), m1, m2, k1, k2, c1, c2, H, L, v_speed);
end

fprintf('Det maximala stabila tidssteget Delta_t_max är: %.5f\n', Delta_t_max);

% Plotta resultat
figure;
plot(t_span, v(1, :), 'DisplayName', 'z1'); % Antag att v(1, :) är z1
hold on;
plot(t_span, v(2, :), 'DisplayName', 'z2'); % Antag att v(2, :) är z2
hold off;
title('Styvare * 100 med Euler framåt');
xlabel('Tid (s)');
ylabel('Förskjutning (m)');
legend;






