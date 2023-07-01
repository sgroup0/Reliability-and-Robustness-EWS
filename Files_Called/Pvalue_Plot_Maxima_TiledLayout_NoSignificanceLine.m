%% PLOT EWS MEASURES
marker_size = 8;
line_width = 1.5;
% Var
figure_counter = figure_counter + 1;
figure(figure_counter);
tiledlayout('flow');
nexttile;

% If H == 1 that is positive slope then plot green dot
selection_H_1 = (H_list_var_maxima == 1);
selection_H_not_1 = (H_list_var_maxima ~= 1);
hold on
plot(mu_list(selection_H_1), p_list_var_maxima(selection_H_1), 'DisplayName', 'Pval maxima Var', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Green2, 'MarkerEdgeColor', PS.Green4);
plot(mu_list(selection_H_not_1), p_list_var_maxima(selection_H_not_1), 'LineStyle', 'none', 'LineWidth', line_width, 'Marker', 'x', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Red2, 'MarkerEdgeColor', PS.Red4);

% yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
xline(0.0023, '--k', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.0023');
legend();

% Sk
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;

% If H == 1 that is positive slope then plot green dot
selection_H_1 = (H_list_sk_maxima == 1);
selection_H_not_1 = (H_list_sk_maxima ~= 1);
hold on
plot(mu_list(selection_H_1), p_list_sk_maxima(selection_H_1), 'DisplayName', 'Pval maxima Sk', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Green2, 'MarkerEdgeColor', PS.Green4);
plot(mu_list(selection_H_not_1), p_list_sk_maxima(selection_H_not_1), 'LineStyle', 'none', 'LineWidth', line_width, 'Marker', 'x', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Red2, 'MarkerEdgeColor', PS.Red4);

% yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
xline(0.0023, '--k', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.0023');
legend();

% Kr
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;

% If H == 1 that is positive slope then plot green dot
selection_H_1 = (H_list_kr_maxima == 1);
selection_H_not_1 = (H_list_kr_maxima ~= 1);
hold on
plot(mu_list(selection_H_1), p_list_kr_maxima(selection_H_1), 'DisplayName', 'Pval maxima Kr', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Green2, 'MarkerEdgeColor', PS.Green4);
plot(mu_list(selection_H_not_1), p_list_kr_maxima(selection_H_not_1), 'LineStyle', 'none', 'LineWidth', line_width, 'Marker', 'x', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Red2, 'MarkerEdgeColor', PS.Red4);

% yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
xline(0.0023, '--k', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.0023');
legend();

% AC Lag 1 second
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;

% If H == 1 that is positive slope then plot green dot
selection_H_1 = (H_list_AC_maxima == 1);
selection_H_not_1 = (H_list_AC_maxima ~= 1);
hold on
plot(mu_list(selection_H_1), p_list_AC_maxima(selection_H_1), 'DisplayName', 'Pval maxima AC', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Green2, 'MarkerEdgeColor', PS.Green4);
plot(mu_list(selection_H_not_1), p_list_AC_maxima(selection_H_not_1), 'LineStyle', 'none', 'LineWidth', line_width, 'Marker', 'x', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Red2, 'MarkerEdgeColor', PS.Red4);

% yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
xline(0.0023, '--k', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.0023');
legend();

% Hurst
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;

% If H == 1 that is positive slope then plot green dot
selection_H_1 = (H_list_H_maxima == 1);
selection_H_not_1 = (H_list_H_maxima ~= 1);
hold on
plot(mu_list(selection_H_1), p_list_H_maxima(selection_H_1), 'DisplayName', 'Pval maxima H', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Green2, 'MarkerEdgeColor', PS.Green4);
plot(mu_list(selection_H_not_1), p_list_H_maxima(selection_H_not_1), 'LineStyle', 'none', 'LineWidth', line_width, 'Marker', 'x', 'MarkerSize', marker_size, 'MarkerFaceColor', PS.Red2, 'MarkerEdgeColor', PS.Red4);

% yline(significance_value_tau, '--r', '0.05', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.05');
xline(0.0023, '--k', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', '0.0023');
legend();






