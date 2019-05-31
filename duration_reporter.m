function [ actual_durations unique_actual_durations intended_durations] = duration_reporter( tdt2mat_filename )
% [ actual_durations intended_durations] = duration_reporter( tdt2mat_filename )
%   This function returns two vectors, the first is the complete list of
%   all the trial durations in order of all the trials, the second is a
%   unique list of the presumed intended durations.  The durations were
%   suggested to be lumped together if they occurred within 10ms (Justus).
%   For the purposes of this program the durations are rounded to the
%   nearest ms.  The clustering to within 10 ms is done elsewhere.
%
%   This function assumes that the passed filename is in a "data" subfolder

cmd = ['load data/' tdt2mat_filename ';'];
disp(['executing: ' cmd]);
eval(cmd);
%cmd = ['tdt2mat_data = ' tdt2mat_filename(1:end-4) ';']; % gets rid of '.dat'
%eval(cmd)
S_ON = tdt2mat_data.epocs.S_ON.onset(1:2:end); % times of start in seconds
SOFF = tdt2mat_data.epocs.SOFF.onset;% times of stop stimulus in seconds

disp(['Noting the lengths of S_ON and SOFF are respectively ' num2str(length(S_ON)) ', ' num2str(length(SOFF))]);
if length(SOFF) ~= length(S_ON)
    smallest_dim = min(length(SOFF), length(S_ON));
    disp(['************ readjusting to ' num2str(smallest_dim)])
    S_ON=S_ON(1:smallest_dim);
    SOFF=SOFF(1:smallest_dim);
end

% actual_durations is similar to StimDur but not rounded
actual_durations= (SOFF-S_ON).*1000; % by keeping actual durations but
% converting secs to ms
unique_actual_durations = unique(round(actual_durations));
intended_durations = [0] ; % figure out
end

