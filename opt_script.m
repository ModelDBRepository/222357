% opt_script.m
% This script optimizes the rise and fall est_firing_rates

% note that params are [amplitude slope threshold], example:
params = [0.1000    1.0000   -0.1000];

% do the rising half sigmoid first
f=@(params)opt_sig(params,shifted_theta_rise, shifted_est_contr_firing_rate_rise);

%smoothed_shifted_est_contr_firing_rate_rise = smooth(shifted_est_contr_firing_rate_rise, min(floor(length(shifted_est_contr_firing_rate_rise)/20) 200])); % no need to restrict ave
smoothed_shifted_est_contr_firing_rate_rise = smooth(shifted_est_contr_firing_rate_rise, floor(length(shifted_est_contr_firing_rate_rise)/20));
[amplitude, slope, threshold]=est_half_sigmoid(shifted_theta_rise, smoothed_shifted_est_contr_firing_rate_rise);

% the first_sig_params_guess is a starting point for the optimization (the
% starting point is what is meant by first
first_sig_params_guess=[amplitude 2 threshold];

% do an optimization
%[opt_rise_params,fval] = fminunc(f,[.2 2 -.1]);
[opt_rise_params,fval] = fminunc(f,first_sig_params_guess);

% now for the falling half sigmoid
params = [.1 -2 1.5];
f=@(params)opt_sig(params,shifted_theta_fall, shifted_est_contr_firing_rate_fall);

smoothed_shifted_est_contr_firing_rate_fall = smooth(shifted_est_contr_firing_rate_fall, floor(length(shifted_est_contr_firing_rate_fall)/20));
[amplitude, slope, threshold]=est_half_sigmoid(shifted_theta_fall, smoothed_shifted_est_contr_firing_rate_fall);

first_sig_params_guess=[amplitude -2 threshold];

% do another optimization
%[opt_fall_params,fval] = fminunc(f,[.1 -3 1.5]);
[opt_fall_params,fval] = fminunc(f,first_sig_params_guess);

