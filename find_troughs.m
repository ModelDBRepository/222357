% find_troughs.m 20140325
% derived from the find_peaks
% use:
% load tdt2mat_data_53.mat
% a = tdt2mat_data.streams.BRTH.data;
% % chop off the end of the data because drops to 0 then 0 for a few entries
% c = a(1:end-3137);
% % find the frequency of the samples, fs, something like 24.414 kHz
% fs=tdt2mat_data.streams.BRTH.fs;
% % smooth the breathing trace to 10ms bins
% bin_size = floor(0.010 * fs); % the number of samples in a 10ms bin
% breathing_trace = smooth(c, bin_size);
% [peak_times, peak_indicies] = find_peaks(breathing_trace, fs);
%
% note: it is convenient to define t for x coordinates of plots
% t=0:1/fs:(length(c)-1)/fs;

function [breathmin_times,breathmin_indicies]=find_troughs(breathing_trace, fs, p_indicies)
% the algorithm is to find a threshold crossing of 1e-6 slope and then
% call the peak when the slope goes negative again
threshold = 2e-7;
breathmin_times_index=1; % assume there is at least one peak
deriv = diff(breathing_trace);
p_indicies_index=1;
head_start = floor(0.03*fs); % 0.03 is three hundredths of a second
% alter if needed, 723 is .03s, neccessary to pass jitter at peak brth cycle

while p_indicies_index<=length(p_indicies)
    index=p_indicies(p_indicies_index)+ head_start; % this is a method to move past jitter at the 
    if index<length(deriv)
        while deriv(index)<0 
            index = index + 1; % keep going until starts rising (positive)
            if index==length(deriv)
                break
            end
        end
    end
    if index==length(deriv)
        break
    end
    breathmin_times(breathmin_times_index) = (index-1)/fs; % noted div with Shaina
    breathmin_indicies(breathmin_times_index) = index-1; % the previous pt is closest
    breathmin_times_index = breathmin_times_index + 1;
    p_indicies_index = p_indicies_index +1;
end
end % end of function
