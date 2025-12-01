function [Sh, Eh, Ih, Rh, Am, Sm, Im, u1, u2, u3] = simulate_with_control( ...
    t, dt, N, T, R, H, ...
    Lambdah, rho1, gammah, dh, xia, mua, betah, betam, Kc, muh, lambdah, ...
    w1, W2, W3, W4, ...
    Sh_0, Eh_0, Ih_0, Rh_0, Am_0, Sm_0, Im_0, max_iter)

Sh = zeros(1, N); Eh = zeros(1, N); Ih = zeros(1, N); Rh = zeros(1, N);
Am = zeros(1, N); Sm = zeros(1, N); Im = zeros(1, N);
Sh(1) = Sh_0; Eh(1) = Eh_0; Ih(1) = Ih_0; Rh(1) = Rh_0;
Am(1) = Am_0; Sm(1) = Sm_0; Im(1) = Im_0;

u1 = zeros(1, N);
u2 = zeros(1, N);
u3 = zeros(1, N);

for iter = 1:max_iter

    u1_old = u1; %#ok<NASGU>
    u2_old = u2; %#ok<NASGU>
    u3_old = u3; %#ok<NASGU>

    % Forward sweep
    for i = 1:N-1
        x = [Sh(i); Eh(i); Ih(i); Rh(i); Am(i); Sm(i); Im(i)];
        u = [u1(i); u2(i); u3(i)];
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

        Sh(i+1) = Sh(i) + dt * k1(1);
        Eh(i+1) = Eh(i) + dt * k1(2);
        Ih(i+1) = Ih(i) + dt * k1(3);
        Rh(i+1) = Rh(i) + dt * k1(4);
        Am(i+1) = Am(i) + dt * k1(5);
        Sm(i+1) = Sm(i) + dt * k1(6);
        Im(i+1) = Im(i) + dt * k1(7);
    end

    % Backward sweep
    lambda = zeros(7, N);
    lambda(:, N) = [0; 0; w1; 0; 0; 0; w1]; 
    for i = N-1:-1:1
        x = [Sh(i); Eh(i); Ih(i); Rh(i); Am(i); Sm(i); Im(i)];
        u = [u1(i); u2(i); u3(i)];
        idx = i;
        T_now = min(max(T(idx), 17.58), 34.6);
        R_now = max(R(idx), 0); %#ok<NASGU>
        H_now = max(H(idx), 0);
        alpha = max(0.001, 0.4 * exp(-(T_now - 30)^2 / (2 * 4^2)));
        phim = 0.00386 * T_now * (T_now - 16.58) * sqrt(34.61 - T_now) * (H_now / 80); phim = max(phim, 0);
        quad = max(-0.000828 * T_now^2 + 0.0367 * T_now + 0.522, 1e-6);
        mum = -log(quad);

        Nh_now = max(x(1) + x(2) + x(3) + x(4), 1e-6);
        dlambdadt = zeros(7,1);
        dlambdadt(1) = -(-lambda(1,i)*(muh + (1-u(1))*alpha*betah*x(7)/Nh_now) + ...
                         lambda(2,i)*(1-u(1))*alpha*betah*x(7)/Nh_now);
        dlambdadt(2) = -(lambdah*(lambda(3,i) - lambda(2,i)) - muh*lambda(2,i));
        dlambdadt(3) = -(w1 - lambda(3,i)*(gammah + dh + muh + u(2)) + ...
                         lambda(4,i)*(gammah+u(2)) - lambda(6,i)*(1-u(1))*alpha*betam*x(6)/Nh_now + ...
                         lambda(7,i)*(1-u(1))*alpha*betam*x(6)/Nh_now);
        dlambdadt(4) = -(lambda(1,i)*rho1 - lambda(4,i)*(rho1 + muh));
        dlambdadt(5) = -(lambda(5,i)*(phim*(x(6) + x(7))/Kc - xia - mua) + lambda(6,i)*xia);
        dlambdadt(6) = -(lambda(5,i)*phim*(1 - x(5)/Kc) - lambda(6,i)*((1-u(1))*alpha*betam*x(3)/Nh_now + mum + u(3)) + ...
                         lambda(7,i)*(1-u(1))*alpha*betam*x(3)/Nh_now);
        dlambdadt(7) = -(lambda(5,i)*phim*(1 - x(5)/Kc) + (1-u(1))*alpha*betah*(lambda(2,i)*x(1)/Nh_now - (1-u(1))*alpha*betah*lambda(1,i))*x(1)/Nh_now - ...
                         lambda(7,i)*(mum + u(3)));

        lambda(:,i) = lambda(:,i+1) + dt * dlambdadt;
    end

    % Control update
    for i = 1:N
        x = [Sh(i); Eh(i); Ih(i); Rh(i); Am(i); Sm(i); Im(i)];
        Nh_now = max(x(1) + x(2) + x(3) + x(4), 1e-6);
        T_now = min(max(T(i), 17.58), 34.6);
        alpha = max(0.001, 0.4 * exp(-(T_now - 30)^2 / (2 * 4^2)));
        u1(i) = min(1, max(0, ((lambda(2,i)-lambda(1,i)) * (alpha * betah * x(1) * x(7)/Nh_now) + ...
                           (lambda(7,i)-lambda(6,i)) * (alpha * betam * x(6) * x(3)/Nh_now)) / W2));
        u2(i) = min(1, max(0, (lambda(4,i) - lambda(3,i)) * x(3) / W3));
        u3(i) = min(1, max(0, (lambda(6,i) * x(6) + lambda(7,i) * x(7)) / W4));
    end
end
end
