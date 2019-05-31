function [ add_spikes, add_time, add_stim] = angle_spikecount( delta_theta, number_of_breaths, tce_selected, tce, tce_spikes,start_timestamp)
%calculate the number of spikes in the tce_selected, which could be the stim on to stim off
%1)For every phase(theta) of the breath cycle, spikes, time, and number 
%of stimuli are accumulated. 
% Find the angle associated with each spike and add to the number of spikes
% at the angle which is called add_spikes (a matrix the size of the
% number of angles (number_in_one_breath)
% This function returns the count over the breath cycle.`

%calculate initialize matricies to add up quantities in all breath cycles
number_in_one_breath=length(0:delta_theta:2*pi-delta_theta); % number of angles in each breath
add_spikes=zeros(1, number_in_one_breath);
add_time=zeros(1, number_in_one_breath);
add_stim=zeros(1, number_in_one_breath);

%  to generate the std and ste, calculate firing rate throughout the
%  recording (we use tce arbitrarily to supply the length below),
% however then just select the firing rate's per time where
%  tce_selected is equal to 1
tce_fr=zeros(1, length(tce)); % Generate tce_fr here
tce_spike_index=1;
tce_spike_indicies=find(tce_spikes==1);

%while tce_spike_index<=length(tce_spike_indicies)
    
    %inter_spike_interval=tce(tce_spike_indicies(tce_spike_index))-tce(tce_spike_indicies(tce_spike_index-1));
     % find the inter spike interval time
     % of a spike by iterating over all the spike indicies and then
     % converting those indicies into a time by the time compresssed and
     % expanded vector
     %instantaneous_fr=1/inter_spike_interval;
     %tcr_fr_start_index=tce_spike_indicies(tce_spike_index-1); % index of first spike in the tce times
     %tcr_fr_end_index=tce_spike_indicies(tce_spike_index); % index of first spike in the tce times
     %tce_fr(tcr_fr_start_index:tcr_fr_end_index)=instantaneous_fr;
     %tce_spike_index = tce_spike_index + 1;
%end

% the firing rate is left zero before the first spike and after the last
% spike, so the tce_fr is now setup

% reshape the selected tce times, tce_selected, as matrix where each row is
% a breath, then loop over the rows computing the number of times that a
% column is selected and form the average, std, and stderr for that
% column
% prototype for reshaping: reshape_theta=reshape(theta(1:end-1), number_of_breaths, (length(theta)-1)/number_of_breaths);
tce_selected_per_breath=reshape(tce_selected(1:end-1),  (length(tce_selected)-1)/number_of_breaths, number_of_breaths)';
tce_fr_per_breath=reshape(tce_fr(1:end-1), (length(tce_fr)-1)/number_of_breaths, number_of_breaths)';

col_index=1;  % go over each column (angle in the breath cycle) and form the statistics
fr_ave=zeros(1, number_in_one_breath);
fr_std=zeros(1, number_in_one_breath);
fr_ste=zeros(1, number_in_one_breath);
timestamp_line54=datevec(now-start_timestamp);
disp(['after line 54 ' num2str(timestamp_line54(5)) 'm ' num2str(timestamp_line54(6))]);

while col_index< number_in_one_breath; % number_of_angles doesn't include the last number_in_one_breath
    row_index=1; % go over each breath (each row is a breath's data)
    col_stats=zeros(1,1); % reinitialize
    col_stats_index = 1; % redefined for each breath angle
    while row_index<=number_of_breaths
        if tce_selected_per_breath(row_index, col_index)==1
            col_stats(col_stats_index)=tce_fr_per_breath(row_index, col_index);
            col_stats_index = col_stats_index + 1;
        end
        row_index= row_index+1;
    end
    
    fr_ave(col_index)=mean(col_stats);
    fr_std(col_index)=std(col_stats);
    if ~isempty(col_stats)
        fr_ste(col_index)=std(col_stats)./(sqrt(length(col_stats)));
    end
    col_index = col_index + 1;
end

timestamp_line77=datevec(now-start_timestamp);
disp(['after line 77 ' num2str(timestamp_line77(5)) 'm ' num2str(timestamp_line77(6))]);

% recalculate average by another method
% loop over tce_selected, transfering the stimulation numbers, spikes and time
tce_index=1;
theta=0:delta_theta:2*pi*number_of_breaths;
% save_angles=zeros(1,1);
% save_indicies=zeros(1,1);
% si_index=1;
% sa_index=1; % index for diagnostic save_angles

% some ideas to speed up: if a smaller theta matrix was passed to
% searchclosest_fixed it might speed things up
breath_angle_theta = 0:delta_theta:2*pi-delta_theta; % no need to go bigger for angles less than 2PI
last_breath_angle_index = length(breath_angle_theta);
hist_tally = 0;
tally_saved=0;
% there is a more efficient way to compute these ix's
% Note that tce follows the breath angles so that all the breath angles
% are represented.  Once 

%for tce_index=1:length(tce_selected)
% more efficient to use the while statement to avoid calls to
% searchclosest_fixed when continuous points are being found where
% the last adjacent tce_selected is true (1)
% This initialization prevents the last_tce_index of being true in
% case the first data point was tce_selected(1)==1; :
last_tce_index_selected = -1; % when tce_selected is one more than
ix=0; % in case comes in with tce_selected(1)=1
% the last_tce_index then the search_closest doesn't have to be called
% because in that situation ix = ix + 1
% however when ix==number_in_one_breath then it
% has to be reset to 1 which is where we will compute the
% statistics for 0 and 2 * pi
while tce_index <length(tce_selected)
    if (tce_selected(tce_index)==1)
        if tce_index==last_tce_index_selected + 1
            ix = ix + 1;
            if ix>last_breath_angle_index
                ix=1; % wrap around to start of breath angle 0 if made it to 2*pi
            end
            tally_saved=tally_saved+1;
        else
            ang=mod(theta(tce_index), 2*pi); % 0<=angle <2*pi is the breath cycle angle
            %save_angles(sa_index)=ang;
            %sa_index=sa_index+1;
            [ix, cval]=searchclosest_fixed(breath_angle_theta, ang); % finds the index ix of the angle in the first
            %save_indicies(si_index) = ix;
            %si_index = si_index + 1;
            % breath of theta
        end
        add_spikes(ix) = add_spikes(ix) + tce_spikes(tce_index); % accumulated spikes (when 1) at this angle
        add_time(ix)=add_time(ix)+tce(tce_index+1)-tce(tce_index); % accumulated time at this breath angle
        add_stim(ix)=add_stim(ix) +1; % represents the presence of one more stimulation at this breath angle
        hist_tally=hist_tally+1;
        last_tce_index_selected = tce_index; % when tce_selected is true record last index in case checking next one also is true
        % then the adjacent breath angles can be used as there is no need
        % to refind the ix index
    end
    tce_index=tce_index+1;
end
% disp(['number_in_one_breath = ' num2str(number_in_one_breath) ' compared to length(add_spikes) = ' num2str(length(add_spikes))])
disp(['hist_tally = ' num2str(hist_tally) ', tally_saved = ' num2str(tally_saved)])
timestamp_line103=datevec(now-start_timestamp);
disp(['after line 103 ' num2str(timestamp_line103(5)) 'm ' num2str(timestamp_line103(6))]);

% calculate firing rate for all stim presentations, irrespective of their
% onset times with respect to breath, and glomon or glomoff periods
% sc_angle=add_spikes./(add_time.*add_stim);

end

