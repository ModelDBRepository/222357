% opt_sig.m
% 
% optimize sigmoids to rise or fall firing rates
%
% see documentation for more
% http://www.mathworks.com/help/optim/ug/passing-extra-parameters.html
%
% usage:
% for amplitude slope threshold
% init_params = [.2 2 -.1];
% will try to match
% y_vec = shifted_est_contr_firing_rate_rise;
% over the x values
% x_vec = shifted_theta_rise; %these are the x axis values

function error_value = opt_sig(init_params, x_vec, y_vec)
model_y_vec = sigmoid( init_params, x_vec);
delta_vec = model_y_vec - y_vec;
error_value = sum(delta_vec.*delta_vec);

