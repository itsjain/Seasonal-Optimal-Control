clc;
close all;

% Time grid
t0 = 0;
tf = 365; 
dt = 0.1; 
t = t0:dt:tf;
N = length(t);

% Climate profiles
[T, R, H] = calculate_climate_profiles(t);

% Parameters
Lambdah = 15;
rho1     = 0.05;
gammah   = 0.059;
dh       = 0.001;
xia      = 0.343;
mua      = 0.05;
w1       = 2;
betah    = 0.267;
betam    = 0.4;
Kc       = 10000;
muh      = 0.05;
lambdah  = 0.008;
Nh       = 2900; %#ok<NASGU>

W2 = 1; 
W3 = 15;   
W4 = 5; 

% Initial conditions
Sh_0 = 1500;
Eh_0 = 820; 
Ih_0 = 210;  
Rh_0 = 450;  
Am_0 = 27792; 
Sm_0 = 4430;  
Im_0 = 1452;  

% Baseline (no control)
[Sh_nc, Eh_nc, Ih_nc, Rh_nc, Am_nc, Sm_nc, Im_nc] = simulate_no_control( ...
    t, dt, N, T, R, H, ...
    Lambdah, rho1, gammah, dh, xia, mua, betah, betam, Kc, muh, lambdah, ...
    Sh_0, Eh_0, Ih_0, Rh_0, Am_0, Sm_0, Im_0);

% With optimal control
max_iter = 500;
[Sh, Eh, Ih, Rh, Am, Sm, Im, u1, u2, u3] = simulate_with_control( ...
    t, dt, N, T, R, H, ...
    Lambdah, rho1, gammah, dh, xia, mua, betah, betam, Kc, muh, lambdah, ...
    w1, W2, W3, W4, ...
    Sh_0, Eh_0, Ih_0, Rh_0, Am_0, Sm_0, Im_0, max_iter);

% Plots
plot_results(t, Sh_nc, Eh_nc, Ih_nc, Rh_nc, Am_nc, Sm_nc, Im_nc, ...
                Sh,    Eh,    Ih,    Rh,    Am,    Sm,    Im,    ...
                u1, u2, u3);
