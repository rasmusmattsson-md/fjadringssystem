

% Definiera parametrar
m1 = 465; % massa för chassidelen (kg)
m2 = 55; % massa för hjulet (kg)
k1 = 5350; % fjäderkonstant för chassifjädern (N/m)
k2_ref = 136100; % referensvärde för fjäderkonstant för däckfjädern (N/m)
c1 = 310; % dämpningskonstant för chassifjädern (Ns/m)
c2 = 1250; % dämpningskonstant för däckfjädern (Ns/m)
H = 0.27; % höjd på gupp (m)
L = 1.1; % längd på gupp (m)
v_speed = 63/3.6; % hastighet av fordonet (m/s)
v0 = [0; 0; 0; 0]; % initialvillkor

k2 = 100 * k2_ref; % Uppdaterad fjäderkonstant för däckfjädern

% Skapa systemmatrisen A korrekt
A = [0 0 1 0; 0 0 0 1; -k1/m1 k1/m1 -c1/m1 c1/m1; k1/m2 -(k1+k2)/m2 c1/m2 -(c1+c2)/m2];

% Beräkna egenvärden och bestäm tmax och h steg
eigenvalues = eig(A);
tmax = min((-2*real(eigenvalues))./(abs(eigenvalues).^2));
alpha_values = [1, 10, 100]; % Alpha värden att testa
Delta_t = alpha_values * tmax;

% Förbered figuren för plotting
figure;
hold on;

% Loop över olika alpha värden
for alpha = alpha_values
    % Skala tidssteget
    h = alpha * tmax;

    % Antal tidssteg och tidsvektor
    t_end = 1;
    num_steps = ceil(t_end / h);
    t_span = linspace(0, t_end, num_steps);

    % Initialisera tillståndsvektor
    y = zeros(4, num_steps);
    y(:,1) = [0; 0; 0; 0];

    % Skapa matrisen som används vid varje iteration av Euler bakåt
    I = eye(size(A));
    B = I - h * A; % Justera detta för Euler bakåt

    % Euler bakåt-loop
    for i = 1:(num_steps - 1)
        t = t_span(i);

        % Bestäm h och h' baserat på tiden
        if t <= L/v_speed
            h_height = (H/2)*(1 - cos((2*pi*v_speed*t)/L));
            h_prim = ((H*pi*v_speed)/L)*sin((2*pi*v_speed*t)/L);
        else
            h_height = 0;
            h_prim = 0;
        end

        % Beräkna kraftvektorn F
        F = [0; 0; 0; (k2*h_height)/m2 + (c2*h_prim)/m2];
        
        % Uppdatera tillståndet y med Euler bakåt
        y(:, i+1) = B \ (y(:, i) + h * F);
    end

    % Plotta resultat för aktuellt alpha
    plot(t_span, y(1, :), 'DisplayName', sprintf('z1: Alpha %d', alpha));
    plot(t_span, y(2, :), 'DisplayName', sprintf('z2: Alpha %d', alpha));
end

% Formatera plotten
title('Euler bakåt med olika alpha');
xlabel('Tid (s)');
ylabel('Förskjutning (m)');
legend;
ylim([-0.03 0.3]);
hold off;












