%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created with ♡ by Hampus & Rasmus
%% Available under the MIT-license
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function dvdt = suspension_system(t, v, m1, m2, k1, k2, c1, c2, H, L, v_speed)
    % Funktionen beskriver vägprofilen och beräknar derivatan av systemtillståndet
    if t <= L / v_speed
        h = (H / 2) * (1 - cos((2 * pi * v_speed * t) / L));
        hd = (H * pi * v_speed / L) * sin((2 * pi * v_speed * t) / L);
    else
        h = 0;
        hd = 0;
    end

    A = [0, 0, 1, 0;
         0, 0, 0, 1;
         -k1/m1, k1/m1, -c1/m1, c1/m1;
         k1/m2, -(k1+k2)/m2, c1/m2, -(c1+c2)/m2];

    F = [0; 0; 0; k2/m2*h + c2/m2*hd];

    dvdt = A * v + F;
end
