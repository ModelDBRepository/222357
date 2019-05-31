% sigmoid.m
% usage:
% value = sigmoid([amplitude slope threshold], x_vec)
% value = amplitude ./ ( 1 + exp( -slope.*(x_vec-threshold))

function value = sigmoid(params, x_vec)
amplitude = params(1);
slope = params(2);
threshold = params(3);
value = amplitude ./ ( 1 + exp( -slope.*(x_vec - threshold)) );
end

