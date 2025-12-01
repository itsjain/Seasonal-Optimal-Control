function [T, R, H] = calculate_climate_profiles(tspan)
    avg_temp = [19.0, 22.7, 27.1, 29.9, 30.6, 29.5, 28.1, 27.9, 27.6, 26.3, 23.3, 20.1];
    min_temp = [12.9, 16.5, 21.4, 25.5, 26.9, 27.0, 26.1, 25.8, 25.3, 23.0, 18.4, 14.5];
    max_temp = [25.3, 29.0, 33.3, 35.7, 35.4, 33.1, 31.1, 31.0, 30.8, 30.2, 28.5, 25.9];
    R_mean = 135.87; A_R = 168.41; Phi_R = -1.933;
    H_mean = 73.19;  A_H = 13.80;  Phi_H = -2.2994;
    R = R_mean + A_R * sin((2*pi/365) * tspan + Phi_R);
    H = H_mean + A_H * sin((2*pi/365) * tspan + Phi_H);
    T = zeros(size(tspan));
    for j = 1:length(tspan)
        day = mod(tspan(j), 365); 
        if day == 0, day = 365; end
        m = min(ceil(day / (365 / 12)), 12);
        T_center = avg_temp(m);
        amplitude = (max_temp(m) - min_temp(m)) / 2;
        T(j) = T_center + amplitude * sin((2 * pi / 30.4) * (day - 15));
    end
end
