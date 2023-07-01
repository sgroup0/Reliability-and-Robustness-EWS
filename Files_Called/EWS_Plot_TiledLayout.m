%% PLOT EWS MEASURES

% Original timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
tiledlayout('flow');
nexttile;
plot(time, state_timeseries, 'DisplayName', 'TimeSeries');
xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
legend('Location', 'northwest');
% % RMS
% % figure_counter = figure_counter + 1;
% % figure(figure_counter);
% % tiledlayout('flow');
% nexttile;
% plot(time_vector{n}, rms_omega{n}, 'DisplayName', 'RMS');
% xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
% xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
% legend('Location', 'northwest');
% 
% % Rate RMS
% % figure_counter = figure_counter + 1;
% % figure(figure_counter);
% % tiledlayout('flow');
% nexttile;
% plot(time_vector{n}, rms_omega{n}, 'DisplayName', 'Rate RMS');
% xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
% xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
% legend('Location', 'northwest');

% Var
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(time_window_ends, var_timeseries, 'DisplayName', 'Var');
xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
legend('Location', 'northwest');

% Sk
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(time_window_ends, sk_timeseries, 'DisplayName', 'Sk');
xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
legend('Location', 'northwest');

% Kr
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(time_window_ends, kr_timeseries, 'DisplayName', 'Kr');
xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
legend('Location', 'northwest');

% AC Lag 1 second
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(time_window_ends, AC_timeseries, 'DisplayName', 'AC 1s');
xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
legend('Location', 'northwest');

% Hurst
% figure_counter = figure_counter + 1;
% figure(figure_counter);
% tiledlayout('flow');
nexttile;
plot(time_window_ends, H_timeseries, 'DisplayName', 'H');
xline(time_max_rate, '--r', 'Transition', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Transition');
xline(time_bifn, '--r', 'Bifurcation', 'LineWidth', 2.5, 'FontSize', 12, 'DisplayName', 'Bifurcation');
legend('Location', 'northwest');






