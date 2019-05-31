function [ fr_angle, fr_ave, fr_std, fr_ste ] = angle_fr( delta_theta, number_of_breaths, tce_selected, tce, tce_spikes)
%calculate the evoked tce_selected, whcih could be the stim on to stim off
%1)For every phase(theta) of the breath cycle, spikes, time, and number 
%of stimuli are accumulated. 2)the firing rate is then calculated by 
%the accummulated spikes, time, and number of stimuli.
% This function returns the firing rate averaged over the breath cycle.`
% Example use:
% [ fr_angle, fr_ave, fr_std, fr_ste ] = angle_fr( delta_theta, number_of_breaths, tce_selected, tce, tce_spikes)

%calculate initialize matricies to add up quantities in all breath cycles
number_in_one_breath=length(0:delta_theta:2*pi);
add_spikes=zeros(1, number_in_one_breath);
add_time=zeros(1, number_in_one_breath);
add_stim=zeros(1, number_in_one_breath);


%  to generate the std and ste, calculate firing rate throughout the
%  recording (we use tce arbitrarily to supply the length below),
% however then just select the firing rate's per time where
%  tce_selected is equal to 1
tce_fr=zeros(1, length(tce)); % Generate tce_fr here
tce_spike_index=2;
tce_spike_indicies=find(tce_spikes==1);

while tce_spike_index<=length(tce_spike_indicies)
     inter_spike_interval=tce(tce_spike_indicies(tce_spike_index))-tce(tce_spike_indicies(tce_spike_index-1));
     % find the inter spike interval time
     % of a spike by iterating over all the spike indicies and then
     % converting those indicies into a time by the time compresssed and
     % expanded vector
     instantaneous_fr=1/inter_spike_interval;
     tcr_fr_start_index=tce_spike_indicies(tce_spike_index-1); % index of first spike in the tce times
     tcr_fr_end_index=tce_spike_indicies(tce_spike_index); % index of first spike in the tce times
     tce_fr(tcr_fr_start_index:tcr_fr_end_index)=instantaneous_fr;
     tce_spike_index = tce_spike_index + 1;
end

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

% recalculate average by another method
% loop over tce_selected, transfering the stimulation numbers, spikes and time
tce_index=1;
theta=0:delta_theta:2*pi*number_of_breaths;
% save_angles=zeros(1,1);
% save_indicies=zeros(1,1);
% si_index=1;
% sa_index=1; % index for diagnostic save_angles
while tce_index <length(tce_selected)
    if (tce_selected(tce_index)==1)
        ang=mod(theta(tce_index), 2*pi); % 0<=angle <2*pi is the breath cycle angle
        %save_angles(sa_index)=ang;
        %sa_index=sa_index+1;
        [ix, cval]=searchclosest_fixed(theta, ang); % finds the index ix of the angle in the first
        %save_indicies(si_index) = ix;
        %si_index = si_index + 1;
        % breath of theta
        add_spikes(ix) = add_spikes(ix) + tce_spikes(tce_index); % accumulated spikes (when 1) at this angle
        add_time(ix)=add_time(ix)+tce(tce_index+1)-tce(tce_index); % accumulated time at this breath angle
        add_stim(ix)=add_stim(ix) +1; % represents the presence of one more stimulation at this breath angle
    end
    tce_index=tce_index+1;
end

% calculate firing rate for all stim presentations, irrespective of their
% onset times with respect to breath, and glomon or glomoff periods
fr_angle=add_spikes./(add_time.*add_stim);

end

