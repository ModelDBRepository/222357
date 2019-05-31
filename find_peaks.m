% find_peaks.m
% finds the peaks of the breathing traces as in 
% tst2mat_data_53.mat's tdt2mat_data
% 20140205
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

function [peak_times, peak_indicies]=find_peaks(breathing_trace, fs)
% the algorithm is to find a threshold crossing of 1e-6 slope and then
% call the peak when the slope goes negative again
peak_times_index=1; % assume there is at least one peak
deriv = diff(breathing_trace);
delta_index = floor((1/50)*length(deriv)); % ignore the first and last 1/50th
% of the breathing trace derivative where frequently there is a glitch
% to calculate the threshold based on heuristic statistics guess
threshold  = (1/2)*std(deriv(delta_index:end-delta_index));
% The general algorithm is to wait for the breath derivative to climb above this
% threshold, then once it crosses back to zero, those zeroes are where the peak
% of the breaths are.
disp(['calculated a breath derivative threshold of ' num2str(threshold)])
% old method was to gues that this would work for all cases: threshold = 2e-7; % 1e-6;
i=1;
skip_last_point=0;
while i<length(deriv)
    if deriv(i)>threshold
        index=i;
%         while deriv(index)>0
        while index<length(deriv) && deriv(index)>0
            index = index + 1; % keep going until dip negative
            if index>=length(deriv) % delete if crash 4 lines
                skip_last_point =1;
            end
        end
        i=index; % move i to where can find the next peak
        if ~skip_last_point % delete if crash just this line and end
            peak_times(peak_times_index) = (index-1)/fs; % noted div with Shaina
            peak_indicies(peak_times_index) = index-1; % the previous pt is closest
            peak_times_index = peak_times_index + 1;
        end
    end
    i=i+1;
end

% if the time is within 5 ms of last peek keep looking for
% next peak
periods=diff(peak_times);
period_threshold = 0.005; % check for and ignore breaths shorter then these seconds
small_indicies=find(periods < period_threshold); % find less than period_threshold seconds
                                                % "breaths" which are likely not breaths
if ~isempty(small_indicies)
    new_peak_times(1)=peak_times(1); % the first breath peak is always OK because
    new_peak_indicies(1)=peak_indicies(1); % there is no prior breath to be to close to
    new_peak_time_index=2; % check for breath times after the first one
    old_peak_time_index=2;
    while old_peak_time_index<=length(peak_times)
        if periods(old_peak_time_index-1) < period_threshold
          % skip this short breath
           old_peak_time_index = old_peak_time_index + 1; % skip over
        else
          new_peak_times(new_peak_time_index)=peak_times(old_peak_time_index);
          new_peak_indicies(new_peak_time_index)=peak_indicies(old_peak_time_index);
          old_peak_time_index = old_peak_time_index + 1;
          new_peak_time_index = new_peak_time_index + 1;
        end
    end
    peak_times=new_peak_times;
    peak_indicies = new_peak_indicies;
end

          