% est_half_sigmoid.m
% usage:
% [amplitude, slope, threshold]=est_half_sigmoid(xvalues, yvalues);
% returns estimates for the amplitude slope and threshold
% that is intended to be used in optimizations that find
% a minimum amplitude slope and threshold based on 
% an error function
%

function [amplitude, slope, threshold]=est_half_sigmoid(xvalues, yvalues)

% find a guess for the amplitude
amplitude = max(yvalues);
min_amplitude = min(yvalues);
% find a guess for the threshold
if amplitude < 2 * min_amplitude 
    % in this case the yvalues do not drop below
    % half the amplitude so we will approximate
    % sigmoid as a line
    mid_point = floor(length(xvalues)/2);
    threshold = xvalues(mid_point);
    slope = (yvalues(end) - yvalues(1)) / (xvalues(end) - xvalues(1));
else
    % calculate threshold guess
    % otherwise use method that will find a mid amplitude point(s)
    mid_amplitude = mean([amplitude min_amplitude]);
    error_vec = (yvalues - mid_amplitude).*(yvalues - mid_amplitude);
    min_indicies = find(error_vec==min(error_vec));
    ave_index = floor(mean(min_indicies));
    threshold = xvalues(ave_index);
    % calculate slope guess
    index_scale = floor(length(yvalues)/16); % going back and forth 1/16 gives 1/8 of graphs distance
    index_high = ave_index + index_scale;
    index_low  = ave_index - index_scale;
    % if the threshold drifts to the edge of the data set the slope
    % calculation can ask for points out of range unless the indexes are
    % brought back into range by the below
    if index_low < 0
        index_low=1;
    end
    if index_high > length(yvalues)
        index_high=length(yvalues);
    end
    slope = (yvalues(index_high) - yvalues(index_low))./(xvalues(index_high) - xvalues(index_low));
end

