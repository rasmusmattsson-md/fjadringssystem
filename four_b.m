clearvars
clc


m1 = 465;
m2 = 55;
k1_ref = 5350;
k2_ref = 136100;
c1 = 310;
c2 = 1250;
v = 63/3.6;
H = 0.27;
L = 1.1;
k2 = 100*k2_ref; 
y0 = [0; 0; 0; 0];

% Referenslösningen med ode45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t_ref, y_ref] = ode45(@(t, y) suspension_system(t, y, m1, m2, k1_ref, k2, c1, c2, v, H, L), [0 0.1], y0, options);

% Beräkna felet Eh för varje h
hs = logspace(-4, -1, 10); % Steglängder att testa
errors_euler = zeros(size(hs));
errors_trap = zeros(size(hs));

for i = 1:length(hs)
    h = hs(i);
    % Implementera framåt Euler
    y_euler = y0;
    for n = 1:(0.1/h)
        t = (n-1)*h;
        F = force(t, m2, k2, c2, v, H, L);
        y_euler = y_euler + h * (A(m1, m2, k1_ref, k2, c1, c2) * y_euler + F);
    end
    errors_euler(i) = abs(y_euler(1) - y_ref(end,1));
    
    % Implementera den implicita trapezmetoden
    y_trap = y0;
    B = eye(4) - h/2 * A(m1, m2, k1_ref, k2, c1, c2);
    for n = 1:(0.1/h)
        t = n*h;
        F_current = force(t-h, m2, k2, c2, v, H, L);
        F_next = force(t, m2, k2, c2, v, H, L);
        y_trap = B \ (y_trap + h/2 * (F_current + F_next));
    end
    errors_trap(i) = abs(y_trap(1) - y_ref(end,1));
end

% Plotta resultaten på en loglog-plot
loglog(hs, errors_euler, 'b-o', hs, errors_trap, 'r-*');
xlabel('Steglängd h');
ylabel('Fel Eh');
legend('Framåt Euler', 'Implicit trapez');
title('Felkonvergens');

fprintf('Tidssteg\t\tFelet (Euler)\tFelet (Trapez)\tKonvergensordning (Euler)\n');
fprintf('------------------------------------------------------------------------\n');
for i = 2:length(hs) % Börjar på 2 eftersom vi inte kan beräkna konvergensordning för första steget
    p_euler = log(errors_euler(i)/errors_euler(i-1)) / log(hs(i)/hs(i-1));
    fprintf('%e\t%e\t%e\t%f\n', hs(i), errors_euler(i), errors_trap(i), p_euler);
end


% Hjälpfunktioner
function F = force(t, m2, k2, c2, v, H, L)
    if t <= L/v
        h = (H/2)*(1 - cos((2*pi*v*t)/L));
        h_prim = ((H*pi*v)/L)*sin((2*pi*v*t)/L);
    else
        h = 0;
        h_prim = 0;
    end
    F = [0; 0; 0; (k2*h)/m2 + (c2*h_prim)/m2];
end

function A_matrix = A(m1, m2, k1, k2, c1, c2)
    A_matrix = [0 0 1 0; 
                0 0 0 1; 
                -k1/m1 k1/m1 -c1/m1 c1/m1; 
                k1/m2 -(k1+k2)/m2 c1/m2 -(c1+c2)/m2];
end

function dy = suspension_system(t, y, m1, m2, k1, k2, c1, c2, v, H, L)
    dy = A(m1, m2, k1, k2, c1, c2) * y + force(t, m2, k2, c2, v, H, L);
end



