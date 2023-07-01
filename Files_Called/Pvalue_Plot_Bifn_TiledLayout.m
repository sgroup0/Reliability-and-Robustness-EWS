%% PLOT EWS MEASURES

% Original timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
tiledlayout('flow');
nexttile;
plot(mu_list, p_list_rms_bifn, 'DisplayName', 'Pval bifn RMS');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Var
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_var_bifn, 'DisplayName', 'Pval bifn Var');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Sk
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_sk_bifn, 'DisplayName', 'Pval bifn Sk');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Kr
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_kr_bifn, 'DisplayName', 'Pval bifn Kr');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% AC Lag 1 second
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_AC_bifn, 'DisplayName', 'Pval bifn AC 1s');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();

% Hurst
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(mu_list, p_list_H_bifn, 'DisplayName', 'Pval bifn H');
yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
legend();






