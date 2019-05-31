function [ Z ] = circle_smooth(Y, span)
%circle_smooth Cicular smoothing of data Y by span points
%   This function uses the matlab smooth function to smooth circular data Y
%   by a the "span" number of points.  Since it needs span over 2 points on
%   each side of the data it creates a new matrix that has these extra
%   points on both ends and then calculates the smoothed with smooth.
%   Finally it chops off the extra added points and returns the circularly
%   smoothed middle part.  Note that span is meant to be odd (see matlab documentation
%   on smooth). See the circle_smooth_demo and
%   circle_smooth_diag for "exercising" this function.
%
%   Tom Morse 20141004
length_of_input = length(Y);
if length_of_input<2 || span <3
    Z=Y; % don't do anything to small inputs
    % optionally print an error message here
    disp([' :*** circle_smooth warning ***: input length ' num2str(length(Y)) ' or span length ' num2str(length(span)) ' is rather small'])
else
    % disp([' note that mod(length(Y),2) = ' num2str(mod(length_of_input,2))])
    % assume that Y has one row
    if ~mod(span,2) % if even
        disp(['Note: in circle_smooth an even span = ' num2str(span) ' was reduced to ' num2str(span-1)])
        span = span - 1; % this is what the smooth function would do anyway
    end
    delta_index = floor(span/2);
    if delta_index>=length_of_input
        delta_index=length_of_input-1;
        disp(['Warning: the span ' num2str(span) ' was rather large compared to length of the input vector ' num2str(length_of_input)])
    end
    % disp([' length(Y) = ' num2str(length_of_input) ', delta_index = ' num2str(delta_index)])
    new_Y = [Y(end-delta_index+1:end) Y Y(1:delta_index)];
    smoothed_Y = smooth(new_Y, span);
    % disp(['length(smoothed_Y = ' num2str(length(smoothed_Y)) ])
    Z = smoothed_Y(delta_index+1:end-delta_index);
end
end
