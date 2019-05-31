% evoke_act_half.m

% This program will generate matrices for control breath responses by
% examining ISI stimuli as well as experiment data, enabling us to examine
% evoked neural responses. 
clear all
close all

tdt2mat_data_index=1;

filename='tdt2mat_data_20140814C1RandR1.mat'
filename='tdt2mat_data_119.mat'
cmd=['load data/' filename];
eval(cmd)

analyze_breathing

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


%if you plot (tce) the figure shows y axis=time of each phase
%in indicies. 
%if you plot (theta, tce) the x axis=phase angle in radians

%convert the spike times by binning the data into 1ms bins
%binned_spikes

ts=tdt2mat_data.snips.eNeu.ts;

sortcode=tdt2mat_data.snips.eNeu.sortcode;


tce_spikes=zeros(1,length(theta)); % spike (1) or no spike (0)
        % for each corresponding time compressed or expanded (tce)
ts_index=1;
% loop over the spikes and assign them to the nearest tce time
while ts_index<=length(ts)
    if sortcode(ts_index)==1
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
              [i,cv] = searchclosest_fixed(tce,spike_time);
              tce_spikes(i)=1;
 %        end
    end
    ts_index=ts_index+1;
end


% make a vector of times to be able to sort data by specific periods
tce_exper=zeros(1, length(theta));% set to 1 from stimulus on 
                                  % to stimulus off

tce_resp=zeros(1, length(theta)); % set to 1's from stimulus on to a set 
                                  % response period
response_period=0.050;  %50 ms following the S_ON
response_start_delay = 0; % set to 0.050 %to start looking at 50 ms after
    %the start of the stimulus
% tce_isi=zeros(1, length(theta)); %sets to 1's from stim off to stim on, 
                                    % the ISI period

tce_contr=zeros(1, length(theta)); %sets 1's from 50ms after stim 
                                        %off to stim on

%tce_off=zeros(1, length(theta));%examines a 50 ms off response period

%tce_spont=zeros(1, length(theta));%examines all control recording pre 
                                       %intial stim presentations

%tce_brth_period=zeros(1, length(theta));%set response period for a given 
                                            %breath
brth_period_select=(2*pi)/8; %set to /8 you are looking at 1/4 of the brth 
                                %cycle centered @ the brth_ang below
                                %Alternatively if you want to look at %set 
                                %to 25 if you want + or - 25 ms, or 50ms 
                                %window
brth_ang_select=pi; %set to pi if you want the trough, set to 2*pi 
                        %if you want the brust

SOFF=tdt2mat_data.epocs.SOFF.onset;

S_ON=tdt2mat_data.epocs.S_ON.onset(1:2:end);

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
tce_exper_start=zeros(1,smallest_length);
tce_exper_end=zeros(1,smallest_length);
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

for s_index= 1 : smallest_length
    %[i_start,cv] = searchclosest_fixed(tce,S_ON(s_index));
    %%%%%%%%%%%%%%%restricts to glomon%%%%%%%%%%%%%
    %[r,c,v]=find(breath_glomon==s_index);%%%%%%
    %if length(v)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [i_start,cv] = searchclosest_fixed(tce,S_ON(s_index)+response_start_delay); 
        
        [i_end,cv] = searchclosest_fixed(tce,SOFF(s_index));
        [i_end_resp,cv] = searchclosest_fixed(tce,S_ON(s_index)+response_period);
        tce_exper_start(s_index) = i_start;
        tce_exper_end(s_index) = i_end; %entire stimulus S_ON (possibly with response start delay)
        tce_resp_start(s_index) = i_start; %narrow window of time, like 15ms following S_ON
        tce_resp_end(s_index) = i_end_resp; %narrow window of time, like 15ms following S_ON
        
        brth_angle=theta(i_start);% gives brth angle for S_ON
        brth_cycle_angle=mod(brth_angle, 2*pi);%give you the radians for the S_ON 
        brth_ang_low=brth_cycle_angle-brth_period_select;
        brth_ang_high=brth_cycle_angle+brth_period_select;
        
        in_interval=0;%set to 1 if you are including burst activity, 
        %which includes stim that are on the borders aka crossing 0 or 2pi
        if brth_ang_low<=brth_ang_select && brth_ang_select<=brth_ang_high
            in_interval=1;
        end
        
        if brth_ang_low<0
            wrapped_brth_ang=mod(brth_ang_low, 2*pi);
            if wrapped_brth_ang<=brth_ang_select
                in_interval=1;
            end
        end
        
        if brth_ang_high>2*pi
            wrapped_brth_ang=mod(brth_ang_high, 2*pi);
            if wrapped_brth_ang>=brth_ang_select
                in_interval=1;
            end
        end
        
%         if in_interval==1
%             tce_brth_period_start(s_index) = i_start;
%             tce_brth_period_end(s_index) = i_end_resp;
%         end
        if s_index<length(S_ON)
            [i_start_next,cv] = searchclosest_fixed(tce,S_ON(s_index+1));
%             tce_isi_start(s_index) = i_end;
%             tce_isi_end(s_index) = i_start_next;
            [i_end_later,cv] = searchclosest_fixed(tce,SOFF(s_index)+0.050);
            %alter 0.050 if you do not
            %want to examine a different spontaneous
            %window that is not 150 ms
            tce_contr_start(s_index) = i_end_later;
            tce_contr_end(s_index) = i_start_next;
            
            [i_end_off,cv] = searchclosest_fixed(tce,S_ON(s_index+1)-0.150);
            % Change -0.150 to a different value if you
            % do not want a 50ms off response period, ex. -180
            % would give you a 20ms off response period
            %tce_off_start(s_index) = i_end;
            %tce_off_end(s_index) = i_end_off;
        end   
    %end
    % s_index = s_index + 1; % go to next stimulus
end

% now that the start and end indicies have been found in parallel they
% must be set into the tce_ matricies:

for s_index=1:smallest_length
    i_start = tce_exper_start(s_index);
    i_end = tce_exper_end(s_index); %entire stimulus S_ON (possibly with response start delay)
    tce_exper(i_start:i_end) = 1;

    i_start = tce_resp_start(s_index); %narrow window of time, like 15ms following S_ON
    i_end_resp = tce_resp_end(s_index); %narrow window of time, like 15ms following S_ON
    tce_resp(i_start:i_end_resp) = 1;
    
%     i_start = tce_brth_period_start(s_index);
%     i_end_resp = tce_brth_period_end(s_index);
%     if 0~=i_start
%       tce_brth_period(i_start:i_end_resp) = 1;
%     end
%     i_end = tce_isi_start(s_index);
%     i_start_next = tce_isi_end(s_index);
%     tce_isi(i_end:i_start_next) = 1;
if s_index<smallest_length
    i_end_later = tce_contr_start(s_index);
    i_start_next = tce_contr_end(s_index);
    tce_contr(i_end_later:i_start_next) = 1;
end

    %i_end = tce_off_start(s_index);
    %i_end_off = tce_off_end(s_index);
    %tce_off(i_end:i_end_off) = 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%END PARFOR REPLACEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [i_end,cv] = searchclosest_fixed(tce,S_ON(1));
% tce_spont(1:i_end-1)=1;

%To make a figure displaying the intervals examined, timing diagram use the
%below plot commands
% figure
% hold on
% plot(tce_exper+10)
% plot(tce_resp+8)
% plot(tce_isi+6)
% plot(tce_contr+4)
% plot(tce_spont+2)
% plot(tce_off)
% ax=axis



%calculate the evoked tce_exper, whcih is the stim on to stim off
%1)For every phase(theta) of the breath cycle, spikes, time, and number 
%of stimuli are accumulated. 2)the firing rate is then calculated by 
%the accummulated spikes, time, and number of stimuli.

% % this is a prototype for this  function that will be 
% % angle_fr
% % This function returns the firing rate angles
% 
% % initialize matricies to add up quantities in all breath cycles
% number_in_one_breath=length(0:delta_theta:2*pi);
% add_spikes=zeros(1, number_in_one_breath);
% add_time=zeros(1, number_in_one_breath);
% add_stim=zeros(1, number_in_one_breath);
% 
% % loop over tce_exper, transfering the stimulation numbers, spikes and time
% tce_index=1;
% 
% while tce_index <length(tce_exper)
%     if (tce_exper(tce_index)==1)
%         ang=mod(theta(tce_index), 2*pi); % 0<=angle <2*pi is the breath cycle angle
%         [ix, cval]=searchclosest_fixed(theta, ang); % finds the index ix of the angle in the first
%         % breath of theta
%         add_spikes(ix) = add_spikes(ix) + tce_spikes(tce_index); % accumulated spikes (when 1) at this angle
%         add_time(ix)=add_time(ix)+tce(tce_index+1)-tce(tce_index); % accumulated time at this breath angle
%         add_stim(ix)=add_stim(ix) +1; % represents the presence of one more stimulation at this breath angle
%     end
%     tce_index=tce_index+1;
% end
% 
% % calculate firing rate for all stim presentations, irrespective of their
% % onset times with respect to breath, and glomon or glomoff periods
% fr_angle=zeros(1, number_in_one_breath);
% fr_angle=add_spikes./(add_time.*add_stim);


%%generates a FR ave from the entire stimulus duration presented(unresticted stim onsets) across all
%%breath cycle angles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ fr_angle, fr_ave, fr_std, fr_ste ] = angle_fr( delta_theta, number_of_breaths, tce_exper, tce, tce_spikes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% figure
% plot( 0:delta_theta:2*pi,fr_angle)
% 
% figure
% plot(fr_ave)
% title('fr ave')
% figure
% plot(fr_std)
% title('fr std')
% figure
% plot(fr_ste)
% title('fr ste')

%%examines a FR at specific stim that start in a set range of angles (which
%%we set to 50 ms response period, see above)
% [ fr_angle_brth, fr_ave_brth, fr_std_brth, fr_ste_brth ] = angle_fr( delta_theta, number_of_breaths,tce_brth_period, tce, tce_spikes);
% 
% figure
% 
% %plot(theta(1:length(fr_ave)), fr_ave)
% errorbar(theta(1:length(fr_ave_brth)), fr_ave_brth,fr_ste_brth)
% figure
% plot(theta(1:length(fr_ave_brth)), fr_ave_brth)

% fr_ave_evoked=fr_ave_brth-fr_ave_contr;%dervived from adding up spikes
% fr_angle_evoked=fr_angle_brth-fr_angle_contr;%dervived from adding firing rates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%          UNCOMMENT TO MAKE FR ANGLE FIGURE LINES 315-334 and
%%%%%%%%%          above lines 276-278
%[ fr_angle_contr, fr_ave_contr, fr_std_contr, fr_ste_contr ] = angle_fr( delta_theta, number_of_breaths, tce_contr, tce, tce_spikes);%tce_contr is 
%what is used to select fr regions

%[ fr_angle_resp, fr_ave_resp, fr_std_resp, fr_ste_resp ] = angle_fr( delta_theta, number_of_breaths, tce_resp, tce, tce_spikes);

% figure
% title('The control breath response (FR) compared to the response period of all stimids presented throughout the breath cycle')
% hold on
% %plot(theta(1:length(fr_ave)), fr_ave)
% errorbar(theta(1:length(fr_ave_contr)), fr_ave_contr,fr_ste_contr)
% plot(theta(1:length(fr_ave_contr)),fr_ave_contr, 'c')
% plot(theta(1:length(fr_ave_resp)), fr_ave_resp, 'r')
% 
% fr_ave_evoked=fr_ave_resp-fr_ave_contr;
% figure
% title('The evoked response to all stimids presented throughout the breath cycle, evoked response = Resp FR ave- Control FR ave.')
% hold on
% plot(theta(1:length(fr_ave_evoked)), fr_ave_evoked, 'k')
% plot([0 2*pi],[0 0],'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%GENERATIONG .DAT FILES FOR EVOKED FITNESS FUNCTION%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over stimids
stimids=tdt2mat_data.epocs.S_ON.data(1:2:end);

% start time of stimuli
S_ON=tdt2mat_data.epocs.S_ON.onset(1:2:end);

%time of SOFF of stimuli
length_stimids=min(length(S_ON), length(SOFF));
stimids=stimids(1:length_stimids); % helps when the last S_ON doesn't 
                                   % have an end recorded

% time of spikes
ts=tdt2mat_data.snips.eNeu.ts;

stimid_index=1;



% parfor stimid_index = 1:length(stimids)
%     stim_time = S_ON(stimid_index);
%     response_start=0; % start looking for spikes when the stimulus starts
%     response_end = response_period;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%% we expect to look for spikes up to 50 ms after the%%%%%%%%%%% 
%         %%%%  stimulus start  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % orig    [indicies, values]=find(ts>stim_time+response_start & ts<=stim_time+response_end);
%         [values, indicies]=find(tce>stim_time+response_start & tce<=stim_time+response_end);
%         spike_count(stimid_index)=sum(tce_spikes(indicies));
%         % find the angles in the breath cycle corresponding to the selected
%         % stimulus time so that the average control firing rate can be
%         % calculated and then subtracted from the raw_fr
%         % angles_in_breath_cycle=mod(theta(indicies), 2*pi);
%         % firing_rates_in_breath_cycle=fr_ave_contr(angles_in_breath_cycle);
%         % calculate the number of expected spikes by iterating over the 
%         % time points of the stimulus transformed into the breath cycle 
%         % angle and adding up the product of the time intervals and the 
%         % firing rates.
%         est_num_of_spikes=0;
%         indicies_index=1;
%         while indicies_index < length(indicies)
%             delta_t = tce(indicies(indicies_index)+1) - tce(indicies(indicies_index));
%             angle=mod(theta(indicies(indicies_index)), 2*pi);
%             % find the angle_index corresponding to this angle
%             %theta of the index that search closest responds
%             [i_angle,cv] = searchclosest_fixed(theta,angle);
%             est_num_of_spikes = est_num_of_spikes + delta_t * fr_ave_contr(i_angle);
%             indicies_index = indicies_index+1;
%         end
%         %control_firing_rate=mean(firing_rates_in_breath_cycle);
%         evoked_fr(stimid_index)=spike_count(stimid_index)-est_num_of_spikes;
%         %est_num_of_spike is derived from the contr_reps
%         %also calculate the breath angle the stimulus starts at 
%         % to also write to the dat file
%         [i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index));
%         stim_angle(stimid_index)=mod(theta(i_stimid), 2*pi);
%         [i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index)+0.050);
%         stim_angle_plus_50(stimid_index) = mod(theta(i_stimid), 2*pi);
%         
% end

% cmd=['fid=fopen(''data/sort_breath/' filename '_FR50ms20140814.dat'',''w'');'];
% eval(cmd)
% 
% for stimid_index=1:length(stimids)
%     y = [stimids(stimid_index); evoked_fr(stimid_index); stim_angle(stimid_index); stim_angle_plus_50(stimid_index)];%stim_angle_plus_50 is being used to save an off response window
%     fprintf(fid, '%g %g %g %g\n', y);
% end
% fclose(fid);
%plot([breath_ang-range_ang breath_ang+range_ang], [a a])

[add_spikes_contr, add_time_contr, add_stim_contr] = angle_spikecount( delta_theta, number_of_breaths, tce_contr, tce, tce_spikes);
%tce_contr is what is used to select fr regions

[add_spikes_resp, add_time_resp, add_stim_resp] = angle_spikecount( delta_theta, number_of_breaths, tce_resp, tce, tce_spikes);
% figure
% hold
% title('Spikes for control (black) and stimulation (red)')
% plot(theta(1:length(add_spikes_contr)),add_spikes_contr, 'b')
% plot(theta(1:length(add_spikes_resp)), add_spikes_resp, 'r')
% 
% sc_evoked=(add_spikes_resp./(add_stim_resp.*add_time_resp))-(add_spikes_contr./(add_stim_contr.*add_time_contr));
% figure
% title('The evoked response to all stimids presented throughout the breath cycle, evoked response = Resp SC ave- Control SC ave.')
% hold on
% plot(theta(1:length(sc_evoked)), sc_evoked, 'k')
% plot([0 2*pi],[0 0],'b')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%generate dat file to creat fitness image based on spike counts,%%%%%%%%
%%%%%%%rather than above firing rates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('starting spike count dat file')
stimid_index=1;
est_spikes_index_start=zeros(1,length(stimids));
est_spikes_index_end=zeros(1,length(stimids));
for stimid_index = 1:length(stimids)
    stim_time = S_ON(stimid_index);
    response_start=0; % start looking for spikes when the stimulus starts
    response_end = response_period; % we expect to look for spikes up to 50 ms after the stimulus start
% orig    [indicies, values]=find(ts>stim_time+response_start & ts<=stim_time+response_end);
    [values, indicies]=find(tce>stim_time+response_start & tce<=stim_time+response_end);
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
 
    est_num_of_spikes=0;
    indicies_index=1;
    est_sc=zeros(1,length(stimids));
disp ('starting spike count dat file2')
for stimid_index = 1:length(stimids);
    indicies=est_spikes_index_start(stimid_index):est_spikes_index_end(stimid_index);
    for indicies_index=1:length(indicies)-1
        delta_t = tce(indicies(indicies_index)+1) - tce(indicies(indicies_index));
        angle=mod(theta(indicies(indicies_index)), 2*pi);
        % find the angle_index corresponding to this angle
        %theta of the index that search closest responds
        [i_angle,cv] = searchclosest_fixed(theta,angle);
        est_num_of_spikes = est_num_of_spikes + delta_t * add_spikes_contr(i_angle)/(add_time_contr(i_angle)*add_stim_contr(i_angle));
        
    end
    est_sc(stimid_index)=est_num_of_spikes;
    
end
%control_firing_rate=mean(firing_rates_in_breath_cycle);
disp ('starting spike count dat file3')
for stimid_index=1:length(stimids)
    evoked_sc(stimid_index)=spike_count(stimid_index)-est_sc(stimid_index);   
    % also calculate the breath angle the stimulus starts at to also write
    % to the dat file
    [i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index));
    stim_angle(stimid_index)=mod(theta(i_stimid), 2*pi);
    [i_stimid,cv] = searchclosest_fixed(tce,S_ON(stimid_index)+0.050);
    stim_angle_plus_50(stimid_index) = mod(theta(i_stimid), 2*pi);
end
    


cmd=['fid=fopen(''data/sort_breath/' filename '_SC50ms20140814.dat'',''w'');'];
eval(cmd)


for stimid_index=1:length(stimids)
    y = [stimids(stimid_index); evoked_sc(stimid_index); stim_angle(stimid_index); stim_angle_plus_50(stimid_index)];%stim_angle_plus_50 is being used to save an off response window
    fprintf(fid, '%g %g %g %g\n', y);
end
fclose(fid);
%plot([breath_ang-range_ang breath_ang+range_ang], [a a])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Generate figure of a specific response to a stimulus %%%%%%%%%%
%%%%%%%%%%%%%during a particular phase of the breath cycle%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimid_index=1;
% spike_index=1;
% breath_ang=5; %center of range of angle being examined ,if examining the 
%               %entire breath cycle set to pi
% ang_range=1; %half width of the range of angles, if examining the entire 
%              %breath cycle set to pi
% 
% spike_index=1; % increments as long as spikes are found
% %%%%%%%%looping through all stimids and placing them into the tce matrix%%
% 
% 
% 
% while stimid_index <= length(stimids)
%     stim_time = S_ON(stimid_index);
%     response_start=0; % start looking for spikes when the stimulus starts
%     response_end = 0.200; % instead of using response_period, use total 
%     % stimulation time; % we expect to look for spikes up to 50 ms after 
%     % the stimulus start
%     % orig    [indicies, values]=find(ts>stim_time+response_start & ts<=stim_time+response_end);
%     [values, indicies]=find(tce>stim_time+response_start & tce<=stim_time+response_end);
%     if numel(indicies)==0% this if statement is needed if there is stimulus before 
%         % the first peak, ignore it, as we do not know where it occurs
%         % relative to the breath cycle.
%         continue
%     end
%         
%     % calculate start_angle to see if in range of desired angle range
%     start_angle=mod(theta(indicies(1)), 2*pi);% angle of the start of 
%                                               % the stimulus (onset)
%     
%     in_interval=0; % if the start_angle falls in the range set by 
%                     % breath_ang plus or minus ang_range
%     
%     if start_angle>breath_ang-ang_range
%         if start_angle<breath_ang+ang_range % non wrap around test
%             in_interval=1;
%         end
%     end
%     if breath_ang-ang_range<0
%         %# wraps around from below 0 to 2*pi side
%         if start_angle > 2*pi + breath_ang-ang_range
%             in_interval=1;
%         end
%     end
%     if breath_ang+ang_range>2*pi
%         %# wraps around above 2*pi side back to 0
%         if start_angle < breath_ang+ang_range - 2 *pi
%             in_interval=1;
%         end
%     end
%     if in_interval==1
%         start_time = S_ON(stimid_index);
%         end_time = S_ON(stimid_index)+0.200; % assume 200 ms stimulations
%         [rows, cols,spike_times]=find(tdt2mat_data.snips.eNeu.ts>=start_time & tdt2mat_data.snips.eNeu.ts<=end_time);
%         spike_times_index=1;
%         while spike_times_index<=length(spike_times)
%             stimulus_evoked_spikes(spike_index)=spike_times(spike_times_index)-start_time;
%             spike_index=spike_index+1;
%             spike_times_index = spike_times_index + 1;
%         end
%     end
%     stimid_index=stimid_index+1;
% end
%     
% 
% 
% 
% figure
% title('Response profile for stimids presented during a set time of the breath cycle')
% hold on
% hist(stimulus_evoked_spikes, 200)
