function plot_results(t, Sh_nc, Eh_nc, Ih_nc, Rh_nc, Am_nc, Sm_nc, Im_nc, ...
                          Sh, Eh, Ih, Rh, Am, Sm, Im, u1, u2, u3)

figure;
plot(t, Sh_nc, 'r--', 'LineWidth', 2, 'DisplayName', 'Without Control');
hold on;
plot(t, Sh, 'b', 'LineWidth', 2, 'DisplayName', 'With Control');
xlabel('Time (days)');
ylabel('Susceptible Hosts');
title('Figure 1: Population of Susceptible Hosts');
legend('show');
set(gca, 'LineWidth', 1, 'FontSize', 10, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

figure;
plot(t, Eh_nc, 'r--', 'LineWidth', 2, 'DisplayName', 'Without Control');
hold on;
plot(t, Eh, 'b', 'LineWidth', 2, 'DisplayName', 'With Control');
xlabel('Time (days)');
ylabel('Exposed Hosts');
title('Figure 2: Population of Exposed Hosts');
legend('show');

figure;
plot(t, Ih_nc, 'r--', 'LineWidth', 2, 'DisplayName', 'Without Control');
hold on;
plot(t, Ih, 'b', 'LineWidth', 2, 'DisplayName', 'With Control');
xlabel('Time (days)');
ylabel('Infected Hosts');
title('Figure 3: Population of Infected Hosts');
legend('show');

figure;
plot(t, Eh_nc+Ih_nc, 'r--', 'LineWidth', 2, 'DisplayName', 'Without Control');
hold on;
plot(t, Eh+Ih, 'b', 'LineWidth', 2, 'DisplayName', 'With Control');
xlabel('Time (days)');
ylabel('Cumulative Infection');
legend('show');

figure;
plot(t, Am_nc, 'r--', 'LineWidth', 2, 'DisplayName', 'Without Control');
hold on;
plot(t, Am, 'b', 'LineWidth', 2, 'DisplayName', 'With Control');
xlabel('Time (days)');
ylabel('Aquatic Vectors');
title('Figure 5: Population of Aquatic Vectors');
legend('show');

figure;
plot(t, Sm_nc, 'r--', 'LineWidth', 2, 'DisplayName', 'Without Control');
hold on;
plot(t, Sm, 'b', 'LineWidth', 2, 'DisplayName', 'With Control');
xlabel('Time (days)');
ylabel('Susceptible Vectors');
title('Figure 6: Population of Susceptible Vectors');
legend('show');

figure;
plot(t, Im_nc, 'r--', 'LineWidth', 2, 'DisplayName', 'Without Control');
hold on;
plot(t, Im, 'b', 'LineWidth', 2, 'DisplayName', 'With Control');
xlabel('Time (days)');
ylabel('Infected Vectors');
legend('show');

figure;
plot(t, u1, 'b-', 'LineWidth', 2, 'DisplayName', 'u_1');
hold on;
plot(t, u2, 'g-', 'LineWidth', 2, 'DisplayName', 'u_2');
plot(t, u3, 'r-', 'LineWidth', 2, 'DisplayName', 'u_3');
xlabel('Time (days)');
ylabel('Control Profile (Strategy 7)');
legend('show');

end
