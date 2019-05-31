

% evoke_act_half.m

% This program will generate matrices for control breath responses by
% examining ISI stimuli as well as experiment data, enabling us to examine

close all
start_timestamp=now;
% disp('timestamp start time initialized to 0')
tdt2mat_data_index=1;
response_start_delay=0.05 % start looking for spikes when the stimulus starts %%% Read into response_profile_static_var_dur 
response_window_delay=response_start_delay;

response_period_end=0.1 %The end of the response_period in real time, or the response window + response start delay 
response_window=response_period_end-response_start_delay % duration of the response period. 
selected_duration=.05

%% Control period Settings%%

check_time=0.19%The minimum amount of time between stimulations
contr_rp_delay=0.14;

response_period_end_control=check_time-contr_rp_delay

%control_response_period=check_time-contr_rp_delay;


% new method to use shared_file_pointers.dat to find the tank
fid=fopen('data/shared_file_pointers.dat','r');
line1=fgets(fid);
cmd = ['load ' line1(1:end-1) ';'];
eval(cmd);
disp(['got tank from ' cmd]);
filename = line1(6:end -8);
disp(['set filename to ' filename ]);
analyze_breathing
timestamp_line17=datevec(now-start_timestamp);
disp(['after analyze_breathing ' num2str(timestamp_line17(5)) 'm ' num2str(timestamp_line17(6))]);
ts =     tdt2mat_data.snips.eNeu.ts;
sortcode=tdt2mat_data.snips.eNeu.sortcode;
%time of SOFF of stimuli
SOFF  =  tdt2mat_data.epocs.SOFF.onset(1:end);
% start time of stimuli
S_ON  =  tdt2mat_data.epocs.S_ON.onset(1:2:end-2);

% detect if an SOFF occurs before an S_ON in the start of the recordings
if (SOFF(1)<S_ON(1))
    SOFF = SOFF(2:end); % this will re-align the S_ON and SOFF's together
end

[S_ON, SOFF ] = S_ON_SOFF_processor(S_ON, SOFF);

num_of_stims=min([length(S_ON) length(SOFF)]);
stim_durs=SOFF(1:num_of_stims)-S_ON(1:num_of_stims);
% loop over stimids
stimids=tdt2mat_data.epocs.S_ON.data(1:2:end);
%stimids=[1:length(stimids)];
% time of spikes
ts=tdt2mat_data.snips.eNeu.ts;
% clear tdt2mat_data % keep around if desired to debug
heuristic_range=5000; % once the first quantity of some tce_matrix is
% found it is believed that the next quantity index will be within this
% index distance from the first.  Adjust to larger values as necessary
% for example if you had a cell that spiked less than once every 5 seconds
% the heuristic_range would need to be adjusted to make sure it would
% capture the distance from one spike index to the next (can be determined
% by stopping the debugger to type in the command line:
%   [value, indicies]= find(spike_times==1
%   diff_indicies = diff(indicies);
%   max_diff=max(diff_indicies)
%  The new heuristic_range = max_diff +1;
% or perhaps it would be safere to specifiy
% heuristic_range = 2* max_diff for future use if needed
% to expand the range
% Another place where the range might need to be increased is if there was
% a goofy stimulus that stayed on for longer than 5 sec or had inter-stim
% intervals greater than 5 sec.

%T_ave = (p_times(end)-p_times(1))/((length(p_times))-1);
%periods=diff(p_times);

%T_max = max(periods);
%number_of_angles = floor((T_max*1000)+0.5);% the number of phase angles
%  per breath cylce
%  T_max here will cause the maximum amount of time
%  between phase angles to be 1 ms during the
%  longest breath.  For everyother breath the time
%  between angles is less than 1 ms.
peak_to_trough = breathmin_times-p_times;
trough_to_peak = p_times(2:end)-breathmin_times(1:end-1);
biggest_halfbreath=max([peak_to_trough trough_to_peak]);
number_of_angles = floor((biggest_halfbreath*1000)+0.5);% the number of
%  phase angles per half breath cycle
%  biggest_halfbreath here will cause the maximum
%  amount of time
%  between phase angles to be 1 ms during the
%  longest halfbreath.  For everyother halfbreath
%  the time between angles is less than 1 ms.

delta_theta=(pi)/(number_of_angles-1);%number of angles is per half breath
number_of_breaths=length(p_times)-1;
theta=0:delta_theta:2*pi*number_of_breaths;
%the choice of the number of angles determines that each breath cycle ends
%at exactly at a multiply of 2pi, so that statistics can be done at each
%angle of a breath cycle.

% these variables for just one breath:
breath_angle_theta = 0:delta_theta:2*pi-delta_theta; % no need to go bigger for angles less than 2PI
last_breath_angle_index = length(breath_angle_theta);


%form vectors that correspond to the thetas (each breath cycle)
tce=zeros(1,length(theta));%This will make a vector of zeros where tce=time
index=1;                   %compressed and expanded, time mapped onto
%specific angles
breath=1;
%T=p_times(breath+1)-p_times(breath);


% transfer each breath cycle time periods in time_chunk into tce
while breath <= number_of_breaths
    Thalf1=breathmin_times(breath)-p_times(breath);
    Thalf2=p_times(breath+1)-breathmin_times(breath);
    deltaT1=Thalf1/(number_of_angles-1); %calculate the time interval that
    deltaT2=Thalf2/(number_of_angles-1); %accounts for pts in the breath
    %corresponse to the theta_angles
    time_chunk1=p_times(breath):deltaT1:breathmin_times(breath);
    %peak to trough times set up in a seq. order
    angle_index=1;
    tce_index=(breath-1)* 2*(number_of_angles-1);
    while angle_index <= number_of_angles %assigning times to tce vector
        tce(tce_index+angle_index)=time_chunk1(angle_index);
        angle_index=angle_index+1;
    end
    tce_index=tce_index+(angle_index-2);
    time_chunk2=breathmin_times(breath):deltaT2:p_times(breath+1);
    angle_index=1;
    %tce_index=(breath-1)* number_of_angles;
    while angle_index <= number_of_angles
        tce(tce_index+angle_index)=time_chunk2(angle_index);
        angle_index=angle_index+1;
    end
    breath=breath+1;
end

timestamp_line94=datevec(now-start_timestamp);
disp(['after breathfile tce created ' num2str(timestamp_line94(5)) 'm ' num2str(timestamp_line94(6))]);

%if you plot (tce) the figure shows y axis=time of each phase
%in indicies.
%if you plot (theta, tce) the x axis=phase angle in radians

%convert the spike times by binning the data into 1ms bins
%binned_spikes



tce_spikes=zeros(1,length(theta)); % spike (1) or no spike (0)
% for each corresponding time compressed or expanded (tce)
ts_index=1;
first_time=1;  % call searchclosest_fixed first time and searchclosest_fixed_incr after
% loop over the spikes and assign them to the nearest tce time
while ts_index<=length(ts)
    if sortcode(ts_index)==2%set to 1 for test data 
    %   if sortcode(ts_index)==2
        %if sortcode(ts_index)>=1 && sortcode(ts_index)<=2 %switch to proper sortcode
        spike_time=ts(ts_index);
        %         if (spike_time>=tce(1)) && (spike_time<=tce(end))
        %             % assign a spike when spike time falls within time breathing
        %             % recorded
        %             z=abs(tce-spike_time);
        %             [i]=find(min(z)==z);
        %             tce_spikes(i)=1;
        %if spike_time <5
        %    disp([num2str(spike_time) ' =spike_time'])
        %end
        if first_time
            [i,cv] = searchclosest_fixed(tce,spike_time);
            first_time=0;
        else
            [i, cv]= searchclosest_fixed(tce, spike_time);%searchclosest_fixed_incr(tce, spike_time,i, heuristic_range);
            tce_spikes(i)=1;
        end
        %        end
    end
    ts_index=ts_index+1;
end

timestamp_line129=datevec(now-start_timestamp);
disp(['after spike_time assigned ' num2str(timestamp_line129(5)) 'm ' num2str(timestamp_line129(6))]);

% make a vector of times to be able to sort data by specific periods
% tce_exper=zeros(1, length(theta));% set to 1 from stimulus on
% to stimulus off

tce_resp=zeros(1, length(theta)); % set to 1's from stimulus on to a set
% response period
%response_period=0.10 % increased for arduino from 0.050;  %50 ms following the S_ON
%response_start_delay = 0.0; % set to 0.050 %to start looking at 50 ms after
%the start of the stimulus

% tce_isi=zeros(1, length(theta)); %sets to 1's from stim off to stim on,
% the ISI period

tce_contr=zeros(1, length(theta)); %sets 1's from 50ms after stim
%off to stim on

%tce_off=zeros(1, length(theta));%examines a 50 ms off response period

tce_spont=zeros(1, length(theta));%examines all control recording pre
%intial stim presentations

%tce_brth_period=zeros(1, length(theta));%set response period for a given
%breath
brth_period_select=pi; % (2*pi)/8; %set to /8 you are looking at 1/4 of the brth
%cycle centered @ the brth_ang below
%Alternatively if you want to look at %set
%to 25 if you want + or - 25 ms, or 50ms
%window
brth_ang_select=pi; %set to pi if you want the trough, set to 2*pi
%if you want the brust


%anaylze control (the last 150ms of the ISI) and the experimental (stimulus
%presentation)
s_index=1;
%load data/breath_glomon.dat;

%%%%setting 1's at the time from stim on to stim off for many cases, the
%%%%complete
%%%%stimulus, the stimuli in a breath angle range, the "stimuli" for the
%%%%control which is actually the time before the light stimuli.  All of
%%%%these are used to select time in tce matrix and get passed to a
%%%%variable in functions called tce_selected.
%%%%window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Parfor initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%while s_index<=length(S_ON) & s_index<=length(SOFF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smallest_length=min(length(S_ON), length(SOFF));

% these matricies record the start and end indicies that are
% regions of these matricies (without the start_end) to 1
% the start_index is in the first column and the end_index
% is in the second column
%tce_exper_start=zeros(1,smallest_length);
%tce_exper_end=zeros(1,smallest_length);
tce_resp_start=zeros(1,smallest_length);
tce_resp_end=zeros(1,smallest_length);
%tce_brth_period_start=zeros(smallest_length);
%tce_brth_period_end=zeros(smallest_length);

% tce_isi_start=zeros(smallest_length);
% tce_isi_end=zeros(smallest_length);
tce_contr_start=zeros(1,smallest_length);
tce_contr_end=zeros(1,smallest_length);
%tce_off_start=zeros(smallest_length);
%tce_off_end=zeros(smallest_length);






for s_index= 2 : smallest_length-1
    
    
    stim_time = S_ON(s_index);
    soff_time = SOFF(s_index);
    previous_soff_time = SOFF(s_index-1);
    next_stim_time = S_ON(s_index+1);
    dur_time = soff_time - stim_time;
    
    adequate_intervals=(stim_time-previous_soff_time)>=check_time; %20140113%&& (next_stim_time-soff_time)>= check_time;
    if (selected_duration <= dur_time) && adequate_intervals 
        %[i_start,cv] = searchclosest_fixed(tce,S_ON(s_index));
        %%%%%%%%%%%%%%%restricts to glomon%%%%%%%%%%%%%
        %[r,c,v]=find(breath_glomon==s_index);%%%%%%
        %if length(v)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if s_index==1
        [i_start,cv] = searchclosest_fixed(tce,S_ON(s_index)+response_start_delay);
        %        [i_end,cv] = searchclosest_fixed(tce,SOFF(s_index));
        [i_end_resp,cv] = searchclosest_fixed(tce,S_ON(s_index)+response_period_end);
        %     else
        %         [i_start,cv] = searchclosest_fixed(tce,S_ON(s_index)+response_start_delay);%searchclosest_fixed_incr(tce,S_ON(s_index)+response_start_delay, i_start, heuristic_range);
        %
        %         %        [i_end,cv] = searchclosest_fixed(tce,SOFF(s_index));
        %         [i_end_resp,cv] = searchclosest_fixed(tce,S_ON(s_index)+response_period_end); % assumes response_period_end>response_start_delay
        %  %searchclosest_fixed_incr(tce,S_ON(s_index)+response_period_end, i_start, heuristic_range)
        %     end
        
        %        tce_exper_start(s_index) = i_start;
        %        tce_exper_end(s_index) = i_end; %entire stimulus S_ON (possibly with response start delay)
        tce_resp_start(s_index) = i_start; %narrow window of time, like 15ms following S_ON
        tce_resp_end(s_index) = i_end_resp; %narrow window of time, like 15ms following S_ON
        
        if s_index > 1
            %            [i_start_next,cv] = searchclosest_fixed(tce,S_ON(s_index+1));
            %             tce_isi_start(s_index) = i_end;
            %             tce_isi_end(s_index) = i_start_next;
            [i_contr_end,cv] = searchclosest_fixed(tce,SOFF(s_index)+contr_rp_delay);
            %[i_end_later,cv] =
            %searchclosest_fixed(tce,S_ON(s_index+1));original 20150303
            [i_end_later,cv] = searchclosest_fixed(tce,SOFF(s_index)+contr_rp_delay+response_period_end_control);% new 20150303
            %window that is not 150 ms
            tce_contr_start(s_index) = i_contr_end;%i_end_later;
            tce_contr_end(s_index) = i_end_later;
            
            %            [i_end_off,cv] = searchclosest_fixed(tce,S_ON(s_index+1)-0.150);
            % Change -0.150 to a different value if you
            % do not want a 50ms off response period, ex. -180
            % would give you a 20ms off response period
            %tce_off_start(s_index) = i_end;
            %tce_off_end(s_index) = i_end_off;
        end
    end
    % s_index = s_index + 1; % go to next stimulus
end

timestamp_line268=datevec(now-start_timestamp);
disp(['after line268 ' num2str(timestamp_line268(5)) 'm ' num2str(timestamp_line268(6))]);

% now that the start and end indicies have been found in parallel they
% must be set into the tce_ matricies:

for s_index=2:smallest_length-1
    %     i_start = tce_exper_start(s_index);
    %     i_end = tce_exper_end(s_index); %entire stimulus S_ON (possibly with response start delay)
    %     tce_exper(i_start:i_end) = 1;
    
    i_start = tce_resp_start(s_index); %narrow window of time, like 15ms following S_ON
    i_end_resp = tce_resp_end(s_index); %narrow window of time, like 15ms following S_ON
    if i_start*i_end_resp>0 
        tce_resp(i_start:i_end_resp) = 1;
    end
    %     i_start = tce_brth_period_start(s_index);
    %     i_end_resp = tce_brth_period_end(s_index);
    %     if 0~=i_start
    %       tce_brth_period(i_start:i_end_resp) = 1;
    %     end
    %     i_end = tce_isi_start(s_index);
    %     i_start_next = tce_isi_end(s_index);
    %     tce_isi(i_end:i_start_next) = 1;
    if s_index<smallest_length-1
        i_end_later = tce_contr_start(s_index);
        i_start_next = tce_contr_end(s_index);
        if i_end_later*i_start_next>0
            tce_contr(i_end_later:i_start_next) = 1;
        end
    end
    
    %i_end = tce_off_start(s_index);
    %i_end_off = tce_off_end(s_index);
    %tce_off(i_end:i_end_off) = 1;
    
end
%%% generates the tce_spont for pre stimulation control histogram
%%% just make it select the time before the first S_ON
i_start=1;
[i_end ,cv] = searchclosest_fixed(tce,S_ON(1));

tce_spont(i_start:i_end)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%GENERATIONG .DAT FILES FOR EVOKED FITNESS FUNCTION%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

length_stimids=min(length(S_ON), length(SOFF));
stimids=stimids(1:length_stimids); % helps when the last S_ON doesn't
% have an end recorded


stimid_index=1;

timestamp_line314=datevec(now-start_timestamp);
disp(['after line 314 ' num2str(timestamp_line314(5)) 'm ' num2str(timestamp_line314(6))]);

[add_spikes_contr, add_time_contr, add_stim_contr] = angle_spikecount( delta_theta, number_of_breaths, tce_contr, tce, tce_spikes,start_timestamp);
%tce_contr is what is used to select fr regions

[add_spikes_resp, add_time_resp, add_stim_resp] = angle_spikecount( delta_theta, number_of_breaths, tce_resp, tce, tce_spikes,start_timestamp);

[add_spikes_spont, add_time_spont, add_stim_spont] = angle_spikecount( delta_theta, number_of_breaths, tce_spont, tce, tce_spikes,start_timestamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%generate dat file to creat fitness image based on spike counts,%%%%%%%%
%%%%%%%rather than above firing rates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timestamp_line327=datevec(now-start_timestamp);
disp(['after line 327 ' num2str(timestamp_line327(5)) 'm ' num2str(timestamp_line327(6))]);

disp ('starting spike count dat file')

stimid_index=1;
est_spikes_index_start=zeros(1,length(stimids));
est_spikes_index_end=zeros(1,length(stimids));
for stimid_index = 1:length(stimids)
    stim_time = S_ON(stimid_index);
    
    % we expect to look for spikes up to 50 ms after the stimulus start
    % orig    [indicies, values]=find(ts>stim_time+response_start & ts<=stim_time+response_end);
    %[values, indicies]=find(tce>stim_time+response_start & tce<=stim_time+response_end);
    % alternate way to form indicies
    [i_start, cv]=searchclosest_fixed(tce,stim_time+response_start_delay);
    [i_end, cv]=searchclosest_fixed(tce,stim_time+response_period_end_control+response_start_delay);% EDITED 20150609SHAINA
    indicies=i_start:i_end;
    spike_count(stimid_index)=sum(tce_spikes(indicies));
    %%% implement searchclosest
    %[i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index));
    
    
    % find the angles in the breath cycle corresponding to the selected
    % stimulus time so that the average control firing rate can be
    % calculated and then subtracted from the raw_fr
    %angles_in_breath_cycle=mod(theta(indicies), 2*pi);
    %firing_rates_in_breath_cycle=fr_ave_contr(angles_in_breath_cycle);
    % calculate the number of expected spikes by iterating over the time
    % points of the stimulus transformed into the breath cycle angle and
    % adding up the product of the time intervals and the firing rates.
    est_spikes_index_start(stimid_index)=indicies(1);
    est_spikes_index_end(stimid_index)=indicies(end);
end

indicies_index=1;
est_sc=zeros(1,length(stimids));

timestamp_line129=datevec(now-start_timestamp);
disp(['at starting dat file2 ' num2str(timestamp_line129(5)) 'm ' num2str(timestamp_line129(6))]);

disp ('starting spike count dat file2')
for stimid_index = 1:length(stimids);
    indicies=est_spikes_index_start(stimid_index):est_spikes_index_end(stimid_index);
    % once we find the start_angle we can cycle through the
    % (stereotyped) breath cycle angles looping back to 0 angle index when
    % needed.
    est_num_of_spikes = 0; % reset to zero each time
    first_time=1;
    for indicies_index=1:length(indicies)-1
        delta_t = tce(indicies(indicies_index)+1) - tce(indicies(indicies_index));
        angle=mod(theta(indicies(indicies_index)), 2*pi);
        % find the angle_index corresponding to this angle
        %theta of the index that search closest responds
        if first_time
            [i_angle,cv] = searchclosest_fixed(breath_angle_theta,angle);
            %[i_angle,cv] = searchclosest_fixed(theta,angle);
            first_time=0;
        else
            i_angle = i_angle + 1;
            if i_angle > last_breath_angle_index
                i_angle=1; % wrap around to start of breath angle 0 if made it to 2*pi
            end
        end
        %est_num_of_spikes = est_num_of_spikes + delta_t * add_spikes_contr(i_angle)/(add_time_contr(i_angle)*add_stim_contr(i_angle));
        est_num_of_spikes = est_num_of_spikes + delta_t * add_spikes_contr(i_angle)/(add_time_contr(i_angle));
        
        n_value=add_stim_contr(i_angle);
        
%         single_control=
%         
%         add_spikes_contr, add_time_contr, add_stim_contr
    end
    
    est_sc(stimid_index)=est_num_of_spikes;
    %est_sd(stimid_index)=
     
end
%control_firing_rate=mean(firing_rates_in_breath_cycle);

timestamp_line129=datevec(now-start_timestamp);
disp(['at starting dat file 3 ' num2str(timestamp_line129(5)) 'm ' num2str(timestamp_line129(6))]);

disp ('starting spike count dat file3')

timestamp_line380=datevec(now-start_timestamp);
disp(['after line 380 ' num2str(timestamp_line380(5)) 'm ' num2str(timestamp_line380(6))]);




for stimid_index=1:length(stimids)
    evoked_sc(stimid_index)=spike_count(stimid_index)-est_sc(stimid_index);
    % also calculate the breath angle the stimulus starts at to also write
    % to the dat file
    [i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index));
    stim_angle(stimid_index)=mod(theta(i_stimid), 2*pi);
    [i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index)+response_start_delay);
    stim_angle_plus_response_start_delay(stimid_index) = mod(theta(i_stimid), 2*pi);
    [i_stimid,cv] = searchclosest_fixed(tce,SOFF(stimid_index));
    soff_angle(stimid_index)=mod(theta(i_stimid), 2*pi);
    [i_stimid,cv] = searchclosest_fixed(tce,SOFF(stimid_index)+contr_rp_delay);
    soff_angle_plus_contr_rp_delay(stimid_index) = mod(theta(i_stimid), 2*pi);
    [i_stimid,cv] = searchclosest_fixed(tce,SOFF(stimid_index)+contr_rp_delay+response_period_end_control);
    soff_angle_plus_response_period_end_control(stimid_index) = mod(theta(i_stimid), 2*pi);
    [i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index)+response_start_delay+.150);
    stim_angle_off_graphing(stimid_index) = mod(theta(i_stimid), 2*pi);
end


cmd=['fid=fopen(''data/sort_breath/' filename '_20150609Epolar.dat'',''w'');'];
eval(cmd)


for stimid_index=1:length(stimids)
    y = [stimids(stimid_index); evoked_sc(stimid_index); ...
        stim_angle(stimid_index); stim_angle_plus_response_start_delay(stimid_index); ...
        soff_angle(stimid_index); soff_angle_plus_contr_rp_delay(stimid_index);...
        soff_angle_plus_response_period_end_control(stimid_index);spike_count(stimid_index)/response_window;...
        contr_rp_delay;selected_duration;response_start_delay;response_period_end; check_time;...
        stim_angle_off_graphing(stimid_index)]%stim_angle_plus_50 is being used to save an off response window
    fprintf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n', y);
end
fclose(fid);
crlf=[char(13) char(10)];
disp(['wrote stimid, evoked spikes, stim angles file:' crlf 'data/sort_breath/' filename '_20150609Dpolar.dat'])
