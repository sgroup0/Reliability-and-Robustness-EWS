%% PLOT EWS MEASURES

% Original timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
tiledlayout('flow');
nexttile;
plot(mu_list, p_list_rms_maxrate, 'DisplayName', 'Pval maxrate RMS');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Var
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_var_maxrate, 'DisplayName', 'Pval maxrate Var');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Sk
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_sk_maxrate, 'DisplayName', 'Pval maxrate Sk');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Kr
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_kr_maxrate, 'DisplayName', 'Pval maxrate Kr');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% AC Lag 1 second
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_AC_maxrate, 'DisplayName', 'Pval maxrate AC 1s');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Hurst
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_H_maxrate, 'DisplayName', 'Pval maxrate H');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();






