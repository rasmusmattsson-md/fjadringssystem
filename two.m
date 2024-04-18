%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created with ♡ by Hampus & Rasmus
%% Available under the MIT-license
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Definiera parametrar
m1 = 465; % massa för chassidelen (kg)
m2 = 55; % massa för hjulet (kg)
k1 = 5350; % fjäderkonstant för chassifjädern (N/m)
k2 = 136100; % fjäderkonstant för däckfjädern (N/m)
c1 = 310; % dämpningskonstant för chassifjädern (Ns/m)
c2 = 1250; % dämpningskonstant för däckfjädern (Ns/m)
H = 0.27; % höjd på gupp (m)
L = 1.1; % längd på gupp (m)
v_speed = 63/3.6; % hastighet av fordonet (m/s)
v0 = [0; 0; 0; 0]; % initialvillkor





% Del a)
options = odeset('RelTol', 1e-6);
[t_ode45, v_ode45] = ode45(@(t, v) suspension_system(t, v, m1, m2, k1, k2, c1, c2, H, L, v_speed), [0 1], v0, options);

% Plotta för del a)
figure;
plot(t_ode45, v_ode45(:, 1), 'DisplayName', 'z1');
hold on;
plot(t_ode45, v_ode45(:, 2), 'DisplayName', 'z2');
hold off;
title('Numerisk lösning för z1 och z2 med ode45');
xlabel('Tid (s)');
ylabel('Förskjutning (m)');
legend;

% Del b)
dt_ode45 = diff(t_ode45);
figure;
plot(t_ode45(1:end-1), dt_ode45, 'DisplayName', 'Tidssteg för ode45');
title('Tidsstegen för ode45');
xlabel('Tid (s)');
ylabel('Tidssteg (s)');
legend;

% Del c)
dt_values = [5e-3, 5e-4];
figure;
for dt_idx = 1:length(dt_values)
    dt = dt_values(dt_idx);
    t_euler = 0:dt:1; % Tidsintervallet från 0 till 1 sekund
    v_euler = zeros(4, length(t_euler));
    v_euler(:, 1) = v0;
    for i = 1:(length(t_euler)-1)
        v_euler(:, i+1) = v_euler(:, i) + dt * suspension_system(t_euler(i), v_euler(:, i), m1, m2, k1, k2, c1, c2, H, L, v_speed);
    end
    plot(t_euler, v_euler(2, :), 'DisplayName', ['Euler Δt = ', num2str(dt)]);
    hold on;
end
plot(t_ode45, v_ode45(:, 2), 'DisplayName', 'ode45', 'LineWidth', 2);
hold off;
title('Jämförelse av z2-lösningar med Euler och ode45');
xlabel('Tid (s)');
ylabel('Förskjutning z2 (m)');
legend;













