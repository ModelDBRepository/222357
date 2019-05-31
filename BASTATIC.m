% Breath_angle_effect.m
% Generate response plots and sort by the breath angle in which the S_ON
% was presented.
%
% This program depends on (reads in) a tank (usually called something like
% tdt2mat_data_something.mat and an evoked spikes breath_info file that
% contains on each line a stimid, number of evoked spikes (experiment -
% control), the start angle that stimulus started at, and an end angle
% (usually the stimulus start time plus 50 ms).
%
% This program also depends on glomon and glomoff dat files that for the
% random block protocol represent the stimid's that illuminated the parent
% glomerulus (glomon) or not (glomoff).  For static files you can either
% prepare a glomon file that contains all the stimids or comment out the
% code that checks for the stimid to be a member of glomon to continue
% adding a stimid to a statistics calculation
RPSTATIC
num_of_breath_increments =12; %24;% 10%40
%response_profile_static_var_dur%This will return the stimids of the stimulaitons that have the proper durations
%select_stimid_adeq_dur
 %Must be the evoked_ program this is the duration of the stimulus;
%control_offset_time=ave_ISI-50;%This is the time post stimulus OFF that contains the off response. This is set so that we do
%not add an off response to the control dataset
control_response_period=response_period;% this is the length of time after a stimulus starts
% that we look for spikes to consider as contributing to a stimulation in
% the "experiment" plot.  We also reuse this number to make the control plot
% from the response_period length of time to look for spikes at starting at
% some other time point (usually SOFF()+0.050, that way the control window
% goes from SOFF()+0.050 to SOFF()+0.050+response_period).
if control_response_period > response_period
    control_period=response_period;
end


%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   Now calculate the triplet histograms and graphs of: 1) control, 2) experiment
%   and 3) evoked = experiment - control
%   The experiment (and control) depend on the response_period = 0.050
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %


%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   First the experiment is computed
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %



%must clear xcenters used in the previous graph
clear xcenters
disp('Calculating experiment')
FR_breath_angles=zeros(1,num_of_breath_increments); % spike counts in breath angles
ste_error_exp = zeros(1, num_of_breath_increments);
exp_num_of_stimuli = zeros(1, num_of_breath_increments);
experiment_stimids={}; % on a per column basis
units_per_trial_exp={};
increment_index = 1;
while increment_index<=num_of_breath_increments
    
    increment_start = 2*pi*(increment_index-1)/num_of_breath_increments;
    increment_end = 2*pi*(increment_index)/num_of_breath_increments;
    xcenters(increment_index)=(increment_start+increment_end)/2-eps; % trying offsetting by eps so 0 is included in next bin
    
    stim_sum=0;
    
    units_per_trial_exp{increment_index}=zeros(1,1);
    
    units_per_trial_index=1;
    per_breath_col_index = 1;
    for index=1:length(breath_info)-1
        % additional check that the stimulation falls under glomon
        if ismember(breath_info(index,1), glomon)% check for glomon under experiment now
            
            stimulation_breath_angle = breath_info(index,exp_breath_info);
            
            if stimulation_breath_angle >= increment_start && stimulation_breath_angle < increment_end
                %         disp([num2str(index) ':found ' num2str(stimulation_breath_angle) ' to be between ' num2str(increment_start) ...
                %             ' and ' num2str(increment_end)])
                
                stimid=breath_info(index,1);
                
                if ismember(stimid, stimids(select_stimid_adeq_dur))
                    
                    %if ismember(stimid, select_breath_angle_stimids)
                    stimid_index=find(stimids==stimid);%find stimid_index in tdt tank to
                    %find the same index in the stimulation time.
                    experiment_stimids{increment_index}(per_breath_col_index) = stimid_index;%n-value that can be used for the ste and std
                    per_breath_col_index = per_breath_col_index + 1;
                    stim_sum=stim_sum+1; %number of stimids per increment
                    
                    stim_time=S_ON(stimid_index)+stim_delay;%now we can use this index to search
                    %through the tdt data tank to retrieve stimulation onset times.
                    
                    %We commonly use to this to sum spikes within response period
                    %'v' in the case is goind to be where we begin to look for
                    %the spike ts stamp
                    %'i' stands for index for x
                    %[i,cv]=searchclosest_fixed(x,v);
                    
                    ap_sum=0;
                    
                    [ts_index,cv]=searchclosest_fixed(ts_raw_sc,stim_time);
                    
                    while ts_index<=length(ts_raw_sc) && ts_raw_sc(ts_index)<stim_time+response_period
                        spike_time = ts_raw_sc(ts_index);% now we need to set an index to loop
                        %over all the spikes within the stimulation window.
                        if spike_time > stim_time
                            
                            %Need to add the proper incrementation of the indices
                            %each_stim_APs(each_stim_APs_index) = spike_time - stim_time;
                            %each_stim_APs_index = each_stim_APs_index + 1;
                            
                            ap_sum=ap_sum+1;
                        end
                        ts_index=ts_index+1;
                    end
                    %calc STE
                    
                    units_per_trial_exp{increment_index}(units_per_trial_index)=ap_sum/response_period;
                    %place this function equal value of each trial
                    units_per_trial_index=units_per_trial_index+1;
                    %need to go through every single trial to generate
                    %the population statistics
                end
            end
        
        end %glomon
        %end
        
    end
    ste_error_exp(increment_index)=std(units_per_trial_exp{increment_index})/sqrt(length(units_per_trial_exp{increment_index}));
    if stim_sum>0
        FR_breath_angles(increment_index)= mean(units_per_trial_exp{increment_index});%increment on the number of bar in the figure
    end
    %     disp(['increment_index = ' num2str(increment_index) ', FR_breath_angles(increment_index) = ' ...
    %         num2str(FR_breath_angles(increment_index))])
    exp_num_of_stimuli(increment_index) = stim_sum;
    increment_index = increment_index + 1;
end

%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
%
%   Secondly the control is computed
%
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %

%%%%%% CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating control')
ste_error_control = zeros(1,num_of_breath_increments);

FR_breath_angles_control = zeros(1,num_of_breath_increments);
control_num_of_stimuli = zeros(1,num_of_breath_increments);
increment_index = 1;
units_per_trial_cont={};
control_stimids={};
while increment_index <= num_of_breath_increments
    increment_start = 2*pi*(increment_index-1)/num_of_breath_increments;
    increment_end =2*pi*(increment_index)/num_of_breath_increments;
    xcenters(increment_index)=(increment_start+increment_end)/2-eps; % offset by eps so 0 is in next bin
    
    stim_sum=0;
    
    units_per_trial_cont{increment_index}=zeros(1,1);
    
    units_per_trial_index=1;
    per_breath_col_index = 1;
    for index=1:length(breath_info)-1
        stimulation_breath_angle = breath_info(index,contr_breath_info); % note index 5 for soff_anlge, index 6 for angle of "SOFF"+ 0.050 secs
        if stimulation_breath_angle >= increment_start && stimulation_breath_angle < increment_end
            
            stimid=breath_info(index,1);
            if ismember(stimid, stimids(select_stimid_adeq_dur))
                stimid_index=find(stimids==stimid);%find stimid_index in tdt tank to
                %find the same index in the stimulation time.
                control_stimids{increment_index}(per_breath_col_index) = stimid_index;
                per_breath_col_index = per_breath_col_index + 1;
                
                %find the same index for the stimulation time's S_ON SOFF.
                stim_sum=stim_sum+1; %number of stimids per increment
                
                
                control_time=SOFF(stimid_index)+control_offset_time; %set here for the control window
                %start time, 50ms is the offset time
                
                %We commonly use to this to sum spikes within response period
                %'v' in the case is goind to be where we begin to look for
                %the spike ts stamp
                %'i' stands for index for x
                %[i,cv]=searchclosest_fixed(x,v);
                
                ap_sum=0;
                
                [ts_index,cv]=searchclosest_fixed(ts_raw_sc,control_time);
                %ts_index<=length(ts_raw_sc) && ts_raw_sc(ts_index)<stim_time+response_period
                while ts_index<=length(ts_raw_sc) && ts_raw_sc(ts_index)<(control_time+control_response_period)%Switch control_response_period SMS 2015 % unless the control response period is
                    % the same as the experiment, the evoked might be difficult
                    % to interpret because very different lengths of times may
                    % be compared
                    spike_time = ts_raw_sc(ts_index);% now we need to set an index to loop
                    %over all the spikes within the stimulation window.
                    if spike_time > control_time
                        
                        %Need to add the proper incrementation of the indices
                        %each_stim_APs(each_stim_APs_index) = spike_time - stim_time;
                        %each_stim_APs_index = each_stim_APs_index + 1;
                        
                        ap_sum=ap_sum+1;
                    end
                    ts_index=ts_index+1;
                end
                %calc STE
                
                units_per_trial_cont{increment_index}(units_per_trial_index)=ap_sum/response_period;
                %place this function equal value of each trial
                units_per_trial_index=units_per_trial_index+1;
                %need to go through every single trial to generate
                %the population statistics
            end
        end
        
    end
    
    ste_error_control(increment_index)=std(units_per_trial_cont{increment_index})/sqrt(length(units_per_trial_cont{increment_index}));
    if stim_sum>0
        FR_breath_angles_control(increment_index)= mean(units_per_trial_cont{increment_index});%increment on the number of bar in the figure
    end
    control_num_of_stimuli(increment_index) = stim_sum;
    increment_index = increment_index + 1;
end


%
% spikes/num_of_stims in each breath_angle bin
%
%
% figure
% hold on
% title('Control Breath Sorted Activity')
% bar(xcenters,FR_breath_angles_control)
% errorbar(xcenters, FR_breath_angles_control, ste_error_control)
% ylabel('spikes/num_of_stims', 'Interpreter', 'none')
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'])
%
% figure
% hold on
% title ('Experimental Breath Sorted Activity')
% bar(xcenters,FR_breath_angles)
% errorbar(xcenters, FR_breath_angles, ste_error_exp)
% ylabel('spikes/num_of_stims', 'Interpreter', 'none')
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'])
%
% figure
% hold on
% title('Evoked Breath Sorted Activity')
% bar(xcenters,(FR_breath_angles-FR_breath_angles_control))
% ylabel('spikes/num_of_stims', 'Interpreter', 'none')
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'])

%
% redone with instantaneous firing rate
%

figure
hold on
title('Control Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
bar(xcenters,FR_breath_angles_control)
errorbar(xcenters, FR_breath_angles_control, ste_error_control)
ylabel('spikes/(num_of_stims*control_response_period) (Hz)', 'Interpreter', 'none','fontsize',label_fontsize)
xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths','fontsize',label_fontsize])

smooth_contr_FR_Breath_angle=circle_smooth(FR_breath_angles_control,3);
[Ycont,Icont] = max(smooth_contr_FR_Breath_angle);
peak_control_angle=(Icont-1)*(2*pi/num_of_breath_increments)+2*pi/(2*num_of_breath_increments)
disp(['Control: peak fr = ' num2str(Ycont) ', at angle = ' num2str(peak_control_angle)]);

FR_breath_angles_control
ste_error_control

figure
hold on
title ('Experimental Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
bar(xcenters,FR_breath_angles)
errorbar(xcenters, FR_breath_angles, ste_error_exp)
ylabel('spikes/(num_of_stims*response_period) (Hz)', 'Interpreter', 'none','fontsize',label_fontsize)
xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'],'fontsize',label_fontsize)

smooth_exp_FR_Breath_angle=circle_smooth(FR_breath_angles,3);
[Yexp,Iexp] = max(smooth_exp_FR_Breath_angle);
peak_exp_angle=(Iexp-1)*(2*pi/num_of_breath_increments)+2*pi/(2*num_of_breath_increments)
disp(['Experimental: peak fr = ' num2str(Yexp) ', at angle = ' num2str(peak_exp_angle)]);
FR_breath_angles
ste_error_exp



% figure
% hold on
% title('Control Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
% bar(xcenters,FR_breath_angles_control./control_response_period)
% errorbar(xcenters, FR_breath_angles_control./control_response_period, ste_error_control./control_response_period)
% ylabel('spikes/(num_of_stims*control_response_period) (Hz)', 'Interpreter', 'none','fontsize',label_fontsize)
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths','fontsize',label_fontsize])


% figure
% hold on
% title ('Experimental Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
% bar(xcenters,FR_breath_angles./response_period)
% errorbar(xcenters, FR_breath_angles./response_period, ste_error_exp./response_period)
% ylabel('spikes/(num_of_stims*response_period) (Hz)', 'Interpreter', 'none','fontsize',label_fontsize)
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'],'fontsize',label_fontsize)

% figure
% hold on
% title('Evoked Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
% bar(xcenters,(FR_breath_angles./response_period - (FR_breath_angles_control./control_response_period)))
% ylabel('spikes/(num_of_stims*response_period) (Hz)', 'Interpreter', 'none','fontsize',label_fontsize)
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'],'fontsize',label_fontsize)

% some diagnostic graphs recording the number of visits to each angle bin

figure
hold on
title('number of control stimuli in breath angle bins','fontsize',label_fontsize)
bar(xcenters, control_num_of_stimuli)
ylabel('num_of_stims (virtual stimulations containing no light)', 'Interpreter', 'none','fontsize',label_fontsize)
xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'],'fontsize',label_fontsize)

figure
hold on
title('number of experiment stimuli in breath angle bins','fontsize',label_fontsize)
bar(xcenters, exp_num_of_stimuli)
ylabel('num_of_stims', 'Interpreter', 'none','fontsize',label_fontsize)
xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths'],'fontsize',label_fontsize)

%Generate print out of p-values comparing the data withing each breath
%angle bin of the control and experimental cells
%Generates P-value for control and experimental breath angle responses to
%understand where a possible phase shift can occured
for index=1:num_of_breath_increments
    [h(index) p(index)]=ttest2(units_per_trial_cont{index},units_per_trial_exp{index});
end
h,p

%


%NOTES
%axis([-0.100    0.3000         0    0.7])
%to see how many stims are analyzes check, size(single_glomon_elements)


% to print the stimids per column in the figure
% for i=1:length(experiment_stimids)
% disp(['Breath col ' num2str(i) ': ' num2str((experiment_stimids{i}))])
% end
% % to print the number of stimids per column in the figure
% for i=1:length(experiment_stimids)
% disp(['Breath col ' num2str(i) ': ' num2str(length(experiment_stimids{i}))])
% end

% Axis Labels
% ax=axis
%
% ax =
%
%      0     7   -50   350
%
% axis([0     7   0   500])

