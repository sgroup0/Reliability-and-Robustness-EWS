%% ESTIMATING KENDALL-TAU VALUES FOR EWS AT DIFFERENT RATES


%% INITIAL SETUP

clear; clc; close all;
PS = PLOT_STANDARDS();
figure_counter = 0;

% Add Files_Called folder to path
addpath('Files_Called');


%% RATES OF CHANGE OF PARAMETER

mu_list = 0.0001: 0.00005: 0.0030;
significance_value_tau = 0.05;


%% ALLOCATE SPACE FOR P-VALUES OF DIFFERENT EWS MEASURES

% Allocate space for timeseries till bifurcation
p_list_rms_bifn = zeros(1, length(mu_list));
p_list_var_bifn = zeros(1, length(mu_list));
p_list_sk_bifn = zeros(1, length(mu_list));
p_list_kr_bifn = zeros(1, length(mu_list));
p_list_AC_bifn = zeros(1, length(mu_list));
p_list_H_bifn = zeros(1, length(mu_list));

% Allocate space for timeseries till max rate of change of rms
p_list_rms_maxrate = zeros(1, length(mu_list));
p_list_var_maxrate = zeros(1, length(mu_list));
p_list_sk_maxrate = zeros(1, length(mu_list));
p_list_kr_maxrate = zeros(1, length(mu_list));
p_list_AC_maxrate = zeros(1, length(mu_list));
p_list_H_maxrate = zeros(1, length(mu_list));

% Allocate space for timeseries till maxima
p_list_var_maxima = zeros(1, length(mu_list));
p_list_sk_maxima = zeros(1, length(mu_list));
p_list_kr_maxima = zeros(1, length(mu_list));
p_list_AC_maxima = zeros(1, length(mu_list));
p_list_H_maxima = zeros(1, length(mu_list));
% Allocate space for H of MMK of timeseries till maxima
H_list_var_maxima = zeros(1, length(mu_list));
H_list_sk_maxima = zeros(1, length(mu_list));
H_list_kr_maxima = zeros(1, length(mu_list));
H_list_AC_maxima = zeros(1, length(mu_list));
H_list_H_maxima = zeros(1, length(mu_list));


%% CALCULATE P-VALUES FOR RATE OF CHANGE OF PARAMETER

% Table of data to be kept track of
Rates = mu_list';
Transient_time = 16 * ones(size(Rates));
Overlap_ratio = (99 / 100) * ones(size(Rates));
Window_size = zeros(size(Rates));
Guess_WS = zeros(size(Rates));
Transition_time = zeros(size(Rates));
Transition_time_transient_removed = zeros(size(Rates));
percentage_WS = zeros(size(Rates));
AC_redline = repmat("W", size(Rates));
H_redline = repmat("W", size(Rates));
TT_mu = zeros(size(Rates));
Guess_TT = zeros(size(Rates));

% Set data for the first Window
% Window_size(1) = 217;
% Transition_time(1) = 738.61;
Window_size(1) = 289;
Transition_time(1) = 738.5;
Transition_time_transient_removed(1) = Transition_time(1) - Transient_time(1);
TT_mu(1) = Transition_time_transient_removed(1) * Rates(1);

%Create table for this data
Check_table = table(Rates, Transient_time, Overlap_ratio, Window_size, Guess_WS, Transition_time, Transition_time_transient_removed, percentage_WS, AC_redline, H_redline, TT_mu, Guess_TT);

for i = 1: length(mu_list)

    fprintf('i = %d\n', i);
    
    % Fix the function inputs
    time_transient = Check_table.Transient_time(i);
    overlap_ratio = Check_table.Overlap_ratio(i);

    if i == 1
        window_size_in_seconds = Check_table.Window_size(1);
    else
        Check_table.Guess_TT(i) = Check_table.TT_mu(i - 1) / Check_table.Rates(i);
        Check_table.Guess_WS(i) = Check_table.Guess_TT(i) * Check_table.percentage_WS(i - 1);
        window_size_in_seconds = Check_table.Guess_WS(i);
    end

    % Call the function
    [transition_time_calc, Ktau_bifn, Ktau_maxrate, Ktau_maxima] = EWS_Ktau_Estimate(i, window_size_in_seconds, time_transient, overlap_ratio);

    % Assign p-value for timeseries till bifurcation
    p_list_rms_bifn(i) = Ktau_bifn.p.rms;
    p_list_var_bifn(i) = Ktau_bifn.p.var;
    p_list_sk_bifn(i) = Ktau_bifn.p.sk;
    p_list_kr_bifn(i) = Ktau_bifn.p.kr;
    p_list_AC_bifn(i) = Ktau_bifn.p.AC;
    p_list_H_bifn(i) = Ktau_bifn.p.H;

    % Assign p-value for timeseries till max rate of change of rms
    p_list_rms_maxrate(i) = Ktau_maxrate.p.rms;
    p_list_var_maxrate(i) = Ktau_maxrate.p.var;
    p_list_sk_maxrate(i) = Ktau_maxrate.p.sk;
    p_list_kr_maxrate(i) = Ktau_maxrate.p.kr;
    p_list_AC_maxrate(i) = Ktau_maxrate.p.AC;
    p_list_H_maxrate(i) = Ktau_maxrate.p.H;

    % Assign p-value for timeseries till maxima of corresponding signal
    p_list_var_maxima(i) = Ktau_maxima.p.var;
    p_list_sk_maxima(i) = Ktau_maxima.p.sk;
    p_list_kr_maxima(i) = Ktau_maxima.p.kr;
    p_list_AC_maxima(i) = Ktau_maxima.p.AC;
    p_list_H_maxima(i) = Ktau_maxima.p.H;
    % Assign H-value for each of the EWS
    H_list_var_maxima(i) = Ktau_maxima.H.var;
    H_list_sk_maxima(i) = Ktau_maxima.H.sk;
    H_list_kr_maxima(i) = Ktau_maxima.H.kr;
    H_list_AC_maxima(i) = Ktau_maxima.H.AC;
    H_list_H_maxima(i) = Ktau_maxima.H.H;
    
    % Give pause to have time to analyze the figures
%     if i < 25
%         i
%         pause(1);
%     else
%         i
%         pause(30);
%     end
%     pause(30);
    

    % Update the data
    Check_table.Transition_time(i) = transition_time_calc;
    Check_table.Transition_time_transient_removed(i) = Check_table.Transition_time(i) - Check_table.Transient_time(i);

    if i ~= 1
        Check_table.Window_size(i) = Check_table.Guess_WS(i);
    end

    Check_table.percentage_WS(i) = Check_table.Window_size(i) / Check_table.Transition_time_transient_removed(i);
    Check_table.TT_mu(i) = Check_table.Transition_time_transient_removed(i) * Check_table.Rates(i);
    
    % Close the generated figures
    close all;
    fprintf('\n');

end

% Assign values to structures to give them as inputs to functions
H_list_maxima.var = H_list_var_maxima;
H_list_maxima.sk = H_list_sk_maxima;
H_list_maxima.kr = H_list_kr_maxima;
H_list_maxima.AC = H_list_AC_maxima;
H_list_maxima.H = H_list_H_maxima;

p_list_maxima.var = p_list_var_maxima;
p_list_maxima.sk = p_list_sk_maxima;
p_list_maxima.kr = p_list_kr_maxima;
p_list_maxima.AC = p_list_AC_maxima;
p_list_maxima.H = p_list_H_maxima;


%% PLOT P-VALUES CHANGE WITH RATE OF CHANGE OF PARAMETER FOR DIFFERENT EWS MEASURES

% Pvalue_Plot_Bifn_TiledLayout
% Pvalue_Plot_Maxrate_TiledLayout

% Print till maxima with significance line
significance_line_print_bool = 1;
figure_counter = Pvalue_Plot_Maxima_TiledLayout(mu_list, H_list_maxima, p_list_maxima, significance_value_tau, significance_line_print_bool, figure_counter);
% Print till maxima without significance line
significance_line_print_bool = 0;
figure_counter = Pvalue_Plot_Maxima_TiledLayout(mu_list, H_list_maxima, p_list_maxima, significance_value_tau, significance_line_print_bool, figure_counter);

% Print till maxima with significance line
significance_line_print_bool = 1;
figure_counter = Pvalue_Plot_Maxima_TiledLayout_OnlyH1(mu_list, H_list_maxima, p_list_maxima, significance_value_tau, significance_line_print_bool, figure_counter);
% Print till maxima without significance line
significance_line_print_bool = 0;
figure_counter = Pvalue_Plot_Maxima_TiledLayout_OnlyH1(mu_list, H_list_maxima, p_list_maxima, significance_value_tau, significance_line_print_bool, figure_counter);


%% SAVE THE FIGURES

% for i = 1: figure_counter
%     figure_name = sprintf('Figures/Pvalue_Plot%d.fig', i);
%     saveas(figure(i), figure_name);
% end





%% CODE FOR CALCULATING KENDALL-TAU AND ASSOCIATED VALUES FUNCTION

function [time_max_rate, Ktau_bifn, Ktau_maxrate, Ktau_maxima] = EWS_Ktau_Estimate(rate_idx, window_size_in_seconds, time_transient, overlap_ratio)
%% INITIAL SETUP

PS = PLOT_STANDARDS();
figure_counter = 0;

%% IMPORT DATA

Pm_bifn = 0.6495;

delta0 = 1;
x0 = cos(delta0);
y0 = sin(delta0);
omega0 = 1.26;
E0 = 1;
Pm0 = .58;

% Time Range details
% nsteps = 1000000;
sampling_rate = 5001;
delta_t = 1 / (sampling_rate - 1);     % the actual formula should be 1 / (sampling_rate - 1), but I use this as an approximation as integer multiple (5000) makes 1 second.
t1 = 0;

Y0 = [x0; y0; omega0; E0; Pm0];

% mu_list = 0.001:0.0005:0.008;
% t2_list = 300 * ones(1, length(mu_list));

% mu_list = 0.0001: 0.00005: 0.0030;
mu_list = 0.0001: 0.00005: 0.0030;

limitcycle_factor = 140 / 100;
Pm_bifn_slope = (0.69 - Pm_bifn) / 0.0023;
Pm_bifn_list = Pm_bifn + Pm_bifn_slope * (mu_list);
t2_list = floor( ((Pm_bifn_list - Pm0) ./ mu_list) * limitcycle_factor );

mu = mu_list(rate_idx);
t2 = t2_list(rate_idx);

time_bifn = (Pm_bifn - Pm0) / mu;

% Load timeseries
filename = sprintf('../Data/Noise5/NoiseOmega5_delta%.2f_omega%.2f_E%.2f_Pm%.4f_mu%.5f_t%.2f_deltaT%.5f_ConstantTimeStep.mat', delta0, omega0, E0, Pm0, mu, t2, delta_t);
Data = load(filename);

tSol = Data.tSol;
YSol = Data.YSol';
xSol = YSol(1, :);
ySol = YSol(2, :);
omegaSol = YSol(3, :);
ESol = YSol(4, :);
PmSol = YSol(5, :);

% Plot omega timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
hold on
plot(tSol, omegaSol);

xlabel('Time');
ylabel('$$\\omega$$', 'Interpreter', 'Latex');

% Set transient time limit
time_transient;

% fprintf('IMPORT DATA\n');
% fprintf('--------------------\n');
% fprintf('t_0 \t\t\t\t\t\t\t= %.2f s\n', t1);
% fprintf('t_f \t\t\t\t\t\t\t= %.2f s\n', t2);
% fprintf('sampling_rate \t\t\t\t\t= %d Hz\n', sampling_rate);
% fprintf('rate_of_parameter_variation \t= %.2f mV/s\n', mu);
% fprintf('parameter_0 \t\t\t\t\t= %.2f\n', PmSol(1));
% fprintf('parameter_f \t\t\t\t\t= %.2f\n', PmSol(end));
% fprintf('time_transient \t\t\t\t\t= %.2f\n', time_transient);
% fprintf('\n\n');


%% CONVERT TO GENERAL VARIABLE NAMES

time = tSol;
state_timeseries = omegaSol;
parameter_variation = PmSol;
sampling_frequency = sampling_rate;
delta_t = delta_t;
parameter_bifurcation = 0.6495;
rate_of_parameter_variation = mu;
bifurcation_time = (parameter_bifurcation - parameter_variation(1)) / rate_of_parameter_variation;


%% REMOVE INTIAL TRANSIENTS AND PLOT NEW TIMESERIES

time_transient = time_transient;
selection_transient = time > time_transient;
time = time(selection_transient);
state_timeseries = state_timeseries(selection_transient);
parameter_variation = parameter_variation(selection_transient);

% Plot transient removed state timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
hold on
plot(time, state_timeseries);


%% LIST OF USEFUL TIMESERIES AND VARIABLES

% fprintf('LIST OF USEFUL TIMESERIES\n');
% fprintf('--------------------\n');
% fprintf('time \t\t\t\t\t\t\t= tSol\n');
% fprintf('state_timeseries \t\t\t\t= omegaSol\n');
% fprintf('parameter_variation \t\t\t= PmSol\n');
% fprintf('sampling_frequency \t\t\t\t= sampling_rate\n');
% fprintf('delta_t \t\t\t\t\t\t= delta_t\n');
% fprintf('parameter_bifurcation \t\t\t= 0.6495\n');
% fprintf('rate_of_parameter_variation \t= mu\n');
% fprintf('\n\n');


%% SET WINDOW DETAILS

overlap_ratio;

window_size_in_seconds;
window_size = floor(window_size_in_seconds / delta_t);

window_step = floor(window_size * (1 - overlap_ratio));
window_step_in_seconds = window_step * delta_t;

% fprintf('SET WINDOW DETAILS\n');
% fprintf('--------------------\n');
% fprintf('window_size \t\t\t\t= %d\n', window_size);
% fprintf('step_size \t\t\t\t\t= %d\n', window_step);
% fprintf('window size in seconds \t\t= %fs\n', window_size_in_seconds);
% fprintf('window step size in seconds \t\t= %fs\n', window_step_in_seconds);
% fprintf('\n\n');


%% GENERATE EWS TIMESERIES

% Generate EWS timeseries for particular window size
[EWS_details, RateRMS_details] = EWS_Timeseries(time, state_timeseries, parameter_variation, delta_t, window_size, window_step);

% Get the timeseries into variables from the structure
time_window_ends = EWS_details.time_window_ends;
parameter_window_ends = EWS_details.parameter_window_ends;
rms_timeseries = EWS_details.rms_timeseries;
var_timeseries = EWS_details.var_timeseries;
sk_timeseries = EWS_details.sk_timeseries;
kr_timeseries = EWS_details.kr_timeseries;
AC_timeseries = EWS_details.AC_timeseries;
H_timeseries = EWS_details.H_timeseries;
H_timeseries_error = EWS_details.H_timeseries_error;

% Return rate of rms details as an object for succinctness
rate_rms_state = RateRMS_details.rate_rms_state;
max_rate_rms_state = RateRMS_details.max_rate_rms_state;
time_vector_rate = RateRMS_details.time_vector_rate;
time_max_rate = RateRMS_details.time_max_rate;
parameter_vector_rate = RateRMS_details.parameter_vector_rate;
parameter_max_rate = RateRMS_details.parameter_max_rate;
rate_transition_time = RateRMS_details.rate_transition_time;


%% TIMESERIES TILL BIFURCATION AND GET NEW TIMESERIES

time_bifn;
selection_bifn_timeseries = time <= time_bifn;
selection_bifn_ews = time_window_ends <= time_bifn;

time_bifn_ts = time(selection_bifn_timeseries);
state_timeseries_bifn = state_timeseries(selection_bifn_timeseries);
parameter_variation_bifn = parameter_variation(selection_bifn_timeseries);

% Get the timeseries into variables from the structure
time_window_ends_bifn = time_window_ends(selection_bifn_ews);
parameter_window_ends_bifn = parameter_window_ends(selection_bifn_ews);
rms_timeseries_bifn = rms_timeseries(selection_bifn_ews);
var_timeseries_bifn = var_timeseries(selection_bifn_ews);
sk_timeseries_bifn = sk_timeseries(selection_bifn_ews);
kr_timeseries_bifn = kr_timeseries(selection_bifn_ews);
AC_timeseries_bifn = AC_timeseries(selection_bifn_ews);
H_timeseries_bifn = H_timeseries(selection_bifn_ews);
H_timeseries_error_bifn = H_timeseries_error(selection_bifn_ews);


%% TIMESERIES TILL MAX RATE AND GET NEW TIMESERIES

time_max_rate;
selection_maxrate_timeseries = time <= time_max_rate;
selection_maxrate_ews = time_window_ends <= time_max_rate;

time_maxrate_ts = time(selection_maxrate_timeseries);
state_timeseries_maxrate = state_timeseries(selection_maxrate_timeseries);
parameter_variation_maxrate = parameter_variation(selection_maxrate_timeseries);

% Get the timeseries into variables from the structure
time_window_ends_maxrate = time_window_ends(selection_maxrate_ews);
parameter_window_ends_maxrate = parameter_window_ends(selection_maxrate_ews);
rms_timeseries_maxrate = rms_timeseries(selection_maxrate_ews);
var_timeseries_maxrate = var_timeseries(selection_maxrate_ews);
sk_timeseries_maxrate = sk_timeseries(selection_maxrate_ews);
kr_timeseries_maxrate = kr_timeseries(selection_maxrate_ews);
AC_timeseries_maxrate = AC_timeseries(selection_maxrate_ews);
H_timeseries_maxrate = H_timeseries(selection_maxrate_ews);
H_timeseries_error_maxrate = H_timeseries_error(selection_maxrate_ews);


%% LIMIT AC AND HURST TILL MAXIMA BEFORE TRANSITION

[~, idx_var_maxima] = max(var_timeseries_maxrate);
[~, idx_sk_maxima] = max(sk_timeseries_maxrate);
[~, idx_kr_maxima] = max(kr_timeseries_maxrate);
[~, idx_AC_maxima] = max(AC_timeseries_maxrate);
[~, idx_H_maxima] = max(H_timeseries_maxrate);

time_var_maxima = time_window_ends_maxrate(idx_var_maxima);
time_sk_maxima = time_window_ends_maxrate(idx_sk_maxima);
time_kr_maxima = time_window_ends_maxrate(idx_kr_maxima);
time_AC_maxima = time_window_ends_maxrate(idx_AC_maxima);
time_H_maxima = time_window_ends_maxrate(idx_H_maxima);

% Var
selection_var_maxima = time_window_ends_maxrate <= time_var_maxima;
time_window_ends_varmaxima = time_window_ends_maxrate(selection_var_maxima);
var_timeseries_maxima = var_timeseries_maxrate(selection_var_maxima);

% Sk
selection_sk_maxima = time_window_ends_maxrate <= time_sk_maxima;
time_window_ends_skmaxima = time_window_ends_maxrate(selection_sk_maxima);
sk_timeseries_maxima = sk_timeseries_maxrate(selection_sk_maxima);

% Kr
selection_kr_maxima = time_window_ends_maxrate <= time_kr_maxima;
time_window_ends_krmaxima = time_window_ends_maxrate(selection_kr_maxima);
kr_timeseries_maxima = kr_timeseries_maxrate(selection_kr_maxima);

% AC
selection_AC_maxima = time_window_ends_maxrate <= time_AC_maxima;
time_window_ends_ACmaxima = time_window_ends_maxrate(selection_AC_maxima);
AC_timeseries_maxima = AC_timeseries_maxrate(selection_AC_maxima);

% H
selection_H_maxima = time_window_ends_maxrate <= time_H_maxima;
time_window_ends_Hmaxima = time_window_ends_maxrate(selection_H_maxima);
H_timeseries_maxima = H_timeseries_maxrate(selection_H_maxima);


%% CALCULATE KENDALL-TAU

% Set significance values
significance_value_tau = 0.05;
significance_value_ac = 0.05;
gpu_shift_critical_size = 520;
print_bool = 0;

default_mmk_assign = 0;

time_ktau_bifn = time_window_ends_bifn;

if length(time_ktau_bifn) > 5
% Calculate Kendall-tau and determine whether to reject or retain null hypothesis
    [tau_bifn.rms, z_bifn.rms, p_bifn.rms, H_bifn.rms] = Modified_MannKendall_test(time_ktau_bifn, rms_timeseries_bifn, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_bifn.var, z_bifn.var, p_bifn.var, H_bifn.var] = Modified_MannKendall_test(time_ktau_bifn, var_timeseries_bifn, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_bifn.sk, z_bifn.sk, p_bifn.sk, H_bifn.sk] = Modified_MannKendall_test(time_ktau_bifn, sk_timeseries_bifn, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_bifn.kr, z_bifn.kr, p_bifn.kr, H_bifn.kr] = Modified_MannKendall_test(time_ktau_bifn, kr_timeseries_bifn, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    
    [tau_bifn.AC, z_bifn.AC, p_bifn.AC, H_bifn.AC] = Modified_MannKendall_test(time_ktau_bifn, AC_timeseries_bifn, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_bifn.H, z_bifn.H, p_bifn.H, H_bifn.H] = Modified_MannKendall_test(time_ktau_bifn, H_timeseries_bifn, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
else
    tau_bifn.rms = default_mmk_assign; z_bifn.rms = default_mmk_assign; p_bifn.rms = default_mmk_assign; H_bifn.rms = default_mmk_assign;
    tau_bifn.var = default_mmk_assign; z_bifn.var = default_mmk_assign; p_bifn.var = default_mmk_assign; H_bifn.var = default_mmk_assign;
    tau_bifn.sk = default_mmk_assign; z_bifn.sk = default_mmk_assign; p_bifn.sk = default_mmk_assign; H_bifn.sk = default_mmk_assign;
    tau_bifn.kr = default_mmk_assign; z_bifn.kr = default_mmk_assign; p_bifn.kr = default_mmk_assign; H_bifn.kr = default_mmk_assign;
    tau_bifn.AC = default_mmk_assign; z_bifn.AC = default_mmk_assign; p_bifn.AC = default_mmk_assign; H_bifn.AC = default_mmk_assign;
    tau_bifn.H = default_mmk_assign; z_bifn.H = default_mmk_assign; p_bifn.H = default_mmk_assign; H_bifn.H = default_mmk_assign;
end

% Calculate Kendall-tau for EWS timeseries till max rate of change of rms
time_ktau_maxrate = time_window_ends_maxrate;

if length(time_ktau_maxrate) > 5
% Calculate Kendall-tau and determine whether to reject or retain null hypothesis
    [tau_maxrate.rms, z_maxrate.rms, p_maxrate.rms, H_maxrate.rms] = Modified_MannKendall_test(time_ktau_maxrate, rms_timeseries_maxrate, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_maxrate.var, z_maxrate.var, p_maxrate.var, H_maxrate.var] = Modified_MannKendall_test(time_ktau_maxrate, var_timeseries_maxrate, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_maxrate.sk, z_maxrate.sk, p_maxrate.sk, H_maxrate.sk] = Modified_MannKendall_test(time_ktau_maxrate, sk_timeseries_maxrate, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_maxrate.kr, z_maxrate.kr, p_maxrate.kr, H_maxrate.kr] = Modified_MannKendall_test(time_ktau_maxrate, kr_timeseries_maxrate, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    
    [tau_maxrate.AC, z_maxrate.AC, p_maxrate.AC, H_maxrate.AC] = Modified_MannKendall_test(time_ktau_maxrate, AC_timeseries_maxrate, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
    [tau_maxrate.H, z_maxrate.H, p_maxrate.H, H_maxrate.H] = Modified_MannKendall_test(time_ktau_maxrate, H_timeseries_maxrate, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
else
    tau_maxrate.rms = default_mmk_assign; z_maxrate.rms = default_mmk_assign; p_maxrate.rms = default_mmk_assign; H_maxrate.rms = default_mmk_assign;
    tau_maxrate.var = default_mmk_assign; z_maxrate.var = default_mmk_assign; p_maxrate.var = default_mmk_assign; H_maxrate.var = default_mmk_assign;
    tau_maxrate.sk = default_mmk_assign; z_maxrate.sk = default_mmk_assign; p_maxrate.sk = default_mmk_assign; H_maxrate.sk = default_mmk_assign;
    tau_maxrate.kr = default_mmk_assign; z_maxrate.kr = default_mmk_assign; p_maxrate.kr = default_mmk_assign; H_maxrate.kr = default_mmk_assign;
    tau_maxrate.AC = default_mmk_assign; z_maxrate.AC = default_mmk_assign; p_maxrate.AC = default_mmk_assign; H_maxrate.AC = default_mmk_assign;
    tau_maxrate.H = default_mmk_assign; z_maxrate.H = default_mmk_assign; p_maxrate.H = default_mmk_assign; H_maxrate.H = default_mmk_assign;
end

% Calculate Kendall-tau for AC and H till their maxima
time_ktau_varmaxima = time_window_ends_varmaxima;
time_ktau_skmaxima = time_window_ends_skmaxima;
time_ktau_krmaxima = time_window_ends_krmaxima;
time_ktau_ACmaxima = time_window_ends_ACmaxima;
time_ktau_Hmaxima = time_window_ends_Hmaxima;

if length(time_ktau_varmaxima) > 5
    [tau_maxima.var, z_maxima.var, p_maxima.var, H_maxima.var] = Modified_MannKendall_test(time_ktau_varmaxima, var_timeseries_maxima, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
else
    tau_maxima.var = default_mmk_assign; z_maxima.var = default_mmk_assign; p_maxima.var = default_mmk_assign; H_maxima.var = default_mmk_assign;
end

if length(time_ktau_skmaxima) > 5
    [tau_maxima.sk, z_maxima.sk, p_maxima.sk, H_maxima.sk] = Modified_MannKendall_test(time_ktau_skmaxima, sk_timeseries_maxima, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
else
    tau_maxima.sk = default_mmk_assign; z_maxima.sk = default_mmk_assign; p_maxima.sk = default_mmk_assign; H_maxima.sk = default_mmk_assign;
end

if length(time_ktau_krmaxima) > 5
    [tau_maxima.kr, z_maxima.kr, p_maxima.kr, H_maxima.kr] = Modified_MannKendall_test(time_ktau_krmaxima, kr_timeseries_maxima, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
else
    tau_maxima.kr = default_mmk_assign; z_maxima.kr = default_mmk_assign; p_maxima.kr = default_mmk_assign; H_maxima.kr = default_mmk_assign;
end

if length(time_ktau_ACmaxima) > 5
    [tau_maxima.AC, z_maxima.AC, p_maxima.AC, H_maxima.AC] = Modified_MannKendall_test(time_ktau_ACmaxima, AC_timeseries_maxima, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
else
    tau_maxima.AC = default_mmk_assign; z_maxima.AC = default_mmk_assign; p_maxima.AC = default_mmk_assign; H_maxima.AC = default_mmk_assign;
end

if length(time_ktau_Hmaxima) > 5
    [tau_maxima.H, z_maxima.H, p_maxima.H, H_maxima.H] = Modified_MannKendall_test(time_ktau_Hmaxima, H_timeseries_maxima, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);
else
    tau_maxima.H = default_mmk_assign; z_maxima.H = default_mmk_assign; p_maxima.H = default_mmk_assign; H_maxima.H = default_mmk_assign;
end

% Save results in structs to return values from function succinctly
Ktau_bifn.tau = tau_bifn;
Ktau_bifn.z = z_bifn;
Ktau_bifn.p = p_bifn;
Ktau_bifn.H = H_bifn;

Ktau_maxrate.tau = tau_maxrate;
Ktau_maxrate.z = z_maxrate;
Ktau_maxrate.p = p_maxrate;
Ktau_maxrate.H = H_maxrate;

Ktau_maxima.tau = tau_maxima;
Ktau_maxima.z = z_maxima;
Ktau_maxima.p = p_maxima;
Ktau_maxima.H = H_maxima;


%% PLOT THE EWS TIMESERIES

EWS_Plot_TiledLayout
EWS_Plot_Bifn_TiledLayout
EWS_Plot_Maxrate_TiledLayout
EWS_Plot_Maxima_TiledLayout


%% SAVE THE FIGURES

% for i = 1: figure_counter
%     figure_name = sprintf('Figures/fig%d.fig', i);
%     saveas(figure(i), figure_name);
% end


end
