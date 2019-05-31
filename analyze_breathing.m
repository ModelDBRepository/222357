% analyze_breathing.

global breathing_trace breathing_trace_derivative breathing_period %period btw two peaks
global breathing_times breathing_indicies %breathing_times are the times of each peak
global tdt2mat_data tdt2mat_data_index
global breathmin_times
%load data/tdt2mat_data_102_103_104.mat; % sets tdt2mat_data, only use if the data file
%has not already been open

a=tdt2mat_data(tdt2mat_data_index).streams.BRTH.data;
c=a;%c=a(1:end-3500); % there are some zeros at the end of the file % for one file we needed 3137
fs=tdt2mat_data(tdt2mat_data_index).streams.BRTH.fs; % Hz, make the sampling rate handy
t=0:1/fs:(length(c)-1)/fs; % a time series to go with breathing trace
num_of_samples = ceil(0.25 * fs); % 0.10 seconds returns about 2440 samples
f=smooth(c,num_of_samples); % 244 samples corresponds to 10 ms at fs sampling rate
breathing_trace=f;
fprime = diff(f); % fprime is proportional to the derivative of breathing rate
%breathing_trace_derivative=fprime;
%figure
%plot(t(1:end-1), fprime)
%title('derivative of breathing - note threshold at 1e-6')

[p_times, p_indicies] = find_peaks(f, fs);
%[p_times, p_indicies] = find_arduino_peaks(f, fs);
% when using find arduino peaks use this _alt instead of regular
% find_troughs

% extra stuff for arduino processing:
% p_times = p_times(1:2:end);
% p_indicies = p_indicies(1:2:end);

%[breathmin_times, breathmin_indicies] = ...
% find_troughs_alt(f, fs, p_indicies); % use for arduino

[breathmin_times, breathmin_indicies] = find_troughs(f, fs, p_indicies);


%hold on
%plot(p_indicies./fs,rectwin(length(p_indicies)).*0,'ro')
%title('breathing trace derivative')

% figure
% plot(t,f)
% hold on
% plot(p_times,f(p_indicies),'rx')
% disp(['length(f)=' num2str(length(f)) ', max(breathmin_indicies) = ' num2str(max(breathmin_indicies))])
usable_breathmin_indicies = find(breathmin_indicies<=length(f));
disp(['length(breathmin_times)=' num2str(length(breathmin_times))])
disp(['length(breathmin_indicies)=' num2str(length(breathmin_indicies))])
disp(['length(usable_breathmin_indicies)=' num2str(length(usable_breathmin_indicies))])
% % use below for incomplete analysis of data from difficult to analyze trace
% plot(breathmin_times(1:length(usable_breathmin_indicies)),f(usable_breathmin_indicies),'bx')
% %plot(breathmin_times,f(breathmin_indicies),'bx')
% %axis([ 0 10 -.002 .003]) % note typing axis and enter on the command line gives current axis size
% title('breathing trace')
breathing_period=diff(p_times);

breathing_times=p_times;
breathing_indicies = p_indicies;
periods=diff(p_times);
breathmin_size=length(breathmin_times);
p_times_size=length(p_times);
if p_times_size ~= breathmin_size % this eliminate the last p-time or breahtmin_time so that
    % there is are equal lengths for both vectors
    min_size=min(breathmin_size, p_times_size);
    breathmin_times=breathmin_times(1:min_size);
    p_times=p_times(1:min_size);
end % there could still be an error if there are two p_times in a row without a breathmin_times between
if length(periods)==length(p_times)
    trough_angle=2*pi*(breathmin_times(1:end-1)-p_times(1:end-1))./periods(1:end-1);
else
    min_end_index = min([length(periods) length(p_times) length(breathmin_times)]);
    trough_angle=2*pi*(breathmin_times(1:min_end_index)-p_times(1:min_end_index))./periods(1:min_end_index);
end
disp(['breathing statistics: mean period = ' num2str(mean(periods)) ', std dev = ' num2str(std(periods))])
disp(['*** look at these --> max period = ' num2str(max(periods)) ', min period = ' num2str(min(periods))])

disp([' trough angle mean = ' num2str(mean(trough_angle)/(2*pi)) '*2*pi, std trough_angle = ' num2str(std(trough_angle)/(2*pi)) '*2*pi'])






