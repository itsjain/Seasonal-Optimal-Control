function [Sh_nc, Eh_nc, Ih_nc, Rh_nc, Am_nc, Sm_nc, Im_nc] = simulate_no_control( ...
    t, dt, N, T, R, H, ...
    Lambdah, rho1, gammah, dh, xia, mua, betah, betam, Kc, muh, lambdah, ...
    Sh_0, Eh_0, Ih_0, Rh_0, Am_0, Sm_0, Im_0)

Sh_nc = zeros(1, N); Eh_nc = zeros(1, N); Ih_nc = zeros(1, N); Rh_nc = zeros(1, N);
Am_nc = zeros(1, N); Sm_nc = zeros(1, N); Im_nc = zeros(1, N);
Sh_nc(1) = Sh_0; Eh_nc(1) = Eh_0; Ih_nc(1) = Ih_0; Rh_nc(1) = Rh_0;
Am_nc(1) = Am_0; Sm_nc(1) = Sm_0; Im_nc(1) = Im_0;

for i = 1:N-1
    x = [Sh_nc(i); Eh_nc(i); Ih_nc(i); Rh_nc(i); Am_nc(i); Sm_nc(i); Im_nc(i)];
    u = [0; 0; 0];
    idx = i;

    T_now = min(max(T(idx), 17.58), 34.6);
    R_now = max(R(idx), 0); %#ok<NASGU>
    H_now = max(H(idx), 0);
    alpha = max(0.001, 0.4 * exp(-(T_now - 30)^2 / (2 * 4^2)));
    phim = 0.00386 * T_now * (T_now - 16.58) * sqrt(34.61 - T_now) * (H_now / 80); phim = max(phim, 0);
    quad = max(-0.000828 * T_now^2 + 0.0367 * T_now + 0.522, 1e-6);
    mum = -log(quad);

    Nh_now = max(x(1) + x(2) + x(3) + x(4), 1e-6);
    k1 = zeros(7,1);
    k1(1) = Lambdah + rho1*x(4) - (1-u(1))*alpha*betah*x(1)*x(7)/Nh_now - muh*x(1);
    k1(2) = (1-u(1))*alpha*betah*x(1)*x(7)/Nh_now - (lambdah + muh)*x(2);
    k1(3) = lambdah*x(2) - (gammah + dh + muh + u(2))*x(3);
    k1(4) = (gammah + u(2))*x(3) - rho1*x(4) - muh*x(4);
    k1(5) = phim*(x(6) + x(7))*(1 - x(5)/Kc) - (xia + mua)*x(5);
    k1(6) = xia*x(5) - (1-u(1))*alpha*betam*x(6)*x(3)/Nh_now - (mum + u(3))*x(6);
    k1(7) = (1-u(1))*alpha*betam*x(6)*x(3)/Nh_now - (mum + u(3))*x(7);

    Sh_nc(i+1) = Sh_nc(i) + dt * k1(1);
    Eh_nc(i+1) = Eh_nc(i) + dt * k1(2);
    Ih_nc(i+1) = Ih_nc(i) + dt * k1(3);
    Rh_nc(i+1) = Rh_nc(i) + dt * k1(4);
    Am_nc(i+1) = Am_nc(i) + dt * k1(5);
    Sm_nc(i+1) = Sm_nc(i) + dt * k1(6);
    Im_nc(i+1) = Im_nc(i) + dt * k1(7);
end
end
