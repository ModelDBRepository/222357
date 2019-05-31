% calc_peak_and_valley_fr_angles.m
% calculate and display the estimated control firing rate and response
% firing rate as a function of breath angle

% this program assumes evoke_act_half_partestdatfile11 or equivalent has
% been run to generate the add_spikes_contr, add_time_contr, add_stim_contr
% and  add_spikes_resp, add_time_resp, add_stim_resp

% est_contr_firing_rate=add_spikes_contr./(add_time_contr.*add_stim_contr);%Shaina
% why this denomenator?
% est_resp_firing_rate=add_spikes_resp./(add_time_resp.*add_stim_resp);

est_contr_firing_rate=add_spikes_contr./(add_time_contr);
est_resp_firing_rate=add_spikes_resp./(add_time_resp);
est_spont_firing_rate=add_spikes_spont./(add_time_spont);

multi_plot = figure;
hold on
subplot(4,2,1) % for 3 rows and two columns plot the 1st graph
hold on
title('unfiltered control firing rate vs breath angle')
plot(breath_angle_theta, est_contr_firing_rate)
ylabel('spikes/(time*num_of_stim) (Hz)', 'Interpreter', 'none')
% [r_s, theta_s]=vector_strength(est_contr_firing_rate,breath_angle_theta);
% plot(theta_s, r_s,'ro')


subplot(4,2,3)
hold on
title('unfiltered response firing rate vs breath angle')
plot(breath_angle_theta, est_resp_firing_rate,'r')
ylabel('spikes/(time*num_of_stim) (Hz)', 'Interpreter', 'none')
% [r_s, theta_s]=vector_strength(est_resp_firing_rate,breath_angle_theta);
% plot(theta_s, r_s,'ro')

subplot(4,2,5)
hold on
title('unfiltered control (blue and response(red) firing rates vs breath angle')
plot(breath_angle_theta, est_resp_firing_rate,'r')
plot(breath_angle_theta, est_contr_firing_rate,'b')
ylabel('spikes/(time*num_of_stim) (Hz)')
subplot(4,2,7)

hold on
title('unfiltered spontaneous (red) firing rates vs breath angle')
plot(breath_angle_theta, est_spont_firing_rate,'r')
ylabel('spikes/(time*num_of_stim) (Hz)')

num_to_ave=floor(.05*length(est_contr_firing_rate)); % was 50; rather than using 50 points choose to use 1/25 of the length
% of the matrix that way it will scale to whatever the breath response is

smoothed_est_contr_firing_rate = circle_smooth(est_contr_firing_rate,num_to_ave);
smoothed_est_resp_firing_rate = circle_smooth(est_resp_firing_rate,num_to_ave);
smoothed_est_spont_firing_rate = circle_smooth(est_spont_firing_rate,num_to_ave);
[Ycont,Icont] = max(smoothed_est_contr_firing_rate);
[Yexp,Iexp] = max(smoothed_est_resp_firing_rate);
[Yspont,Ispont] = max(smoothed_est_spont_firing_rate);

peak_control_angle = breath_angle_theta(Icont);
peak_control_fr = Ycont;
peak_exp_angle = breath_angle_theta(Iexp);
peak_exp_fr = Yexp;
peak_spont_angle = breath_angle_theta(Ispont);
peak_spont_fr = Yspont;
disp(['Control: peak fr = ' num2str(Ycont) ', at angle = ' num2str(breath_angle_theta(Icont))]);
disp(['Experiment: peak fr = ' num2str(Yexp) ', at angle = ' num2str(breath_angle_theta(Iexp))]);
disp(['Spontaneous: peak fr = ' num2str(Yspont) ', at angle = ' num2str(breath_angle_theta(Ispont))]);



subplot(4,2,2)
hold on
title('smoothed control vs breath angle')
plot(breath_angle_theta, smoothed_est_contr_firing_rate)
subplot(4,2,4)
hold on
title('smoothed response firing rate vs breath angle')
plot(breath_angle_theta, smoothed_est_resp_firing_rate, 'r')
subplot(4,2,6)
title('smoothed control (blue and response(red) firing rates vs breath angle')
hold on
plot(breath_angle_theta, smoothed_est_contr_firing_rate)
plot(breath_angle_theta, smoothed_est_resp_firing_rate, 'r')
subplot(4,2,8)
title('smoothed spont (red) firing rates vs breath angle')
hold on
plot(breath_angle_theta, smoothed_est_spont_firing_rate, 'r')

[peak_value, peak_index]=max(smoothed_est_contr_firing_rate);
disp(['peak_value = ' num2str(peak_value) ', angle = ' num2str(2*pi*peak_index/length(smoothed_est_contr_firing_rate))])
num_of_angles_per_cycle=length(est_contr_firing_rate);
for index=1:num_of_angles_per_cycle
    xi(index)=est_contr_firing_rate(index)*cos(2*pi*index/num_of_angles_per_cycle);
    yi(index)=est_contr_firing_rate(index)*sin(2*pi*index/num_of_angles_per_cycle);
end
x_ave=mean(xi);
y_ave=mean(yi);
peak_theta=atan(y_ave/x_ave);

disp(['peak_theta = ' num2str(peak_theta) ' radians = ' num2str(peak_theta/(2*pi)) '*2*pi = ' num2str(peak_theta*360/(2*pi)) ' degrees'])
shifted_theta_axes_start=peak_theta - pi;
shifted_theta_axes_end = peak_theta+pi;
shifted_theta_axes=shifted_theta_axes_start:delta_theta:shifted_theta_axes_end;
shifted_theta_rise=shifted_theta_axes_start:delta_theta:peak_theta;
shifted_theta_fall=peak_theta:delta_theta:shifted_theta_axes_end;

for index=1:length(shifted_theta_rise)
    phi=shifted_theta_rise(index);
    % shift to 0 to 2pi range
    phi = mod(phi, 2*pi);
    % convert to index
    phi_index = floor((phi+delta_theta-1e-9)*num_of_angles_per_cycle/(2*pi));
    shifted_est_contr_firing_rate_rise(index) = est_contr_firing_rate(phi_index);
    post_index_vec(index) = index;
    post_index_vec(index) = phi_index;
end
rise_fig = figure;
plot(shifted_theta_rise, shifted_est_contr_firing_rate_rise)

for index=1:length(shifted_theta_fall)
    phi=shifted_theta_fall(index);
    % shift to 0 to 2pi range
    phi = mod(phi, 2*pi);
    % convert to index
    phi_index = floor((phi+delta_theta-1e-9)*num_of_angles_per_cycle/(2*pi));
    shifted_est_contr_firing_rate_fall(index) = est_contr_firing_rate(phi_index);
    post_index_vec(index) = index;
    post_index_vec(index) = phi_index;
end
fall_fig = figure;
plot(shifted_theta_fall, shifted_est_contr_firing_rate_fall)

opt_script

values= sigmoid(opt_rise_params, shifted_theta_rise);
figure(rise_fig)
hold on
plot(shifted_theta_rise, values, 'r')

values= sigmoid(opt_fall_params, shifted_theta_fall);
figure(fall_fig)
hold on
plot(shifted_theta_fall, values, 'r')

% use
sigmoid_threshold =0.05;
% to find the firing rate angle's where the sigmoid passes through
% this value
% in the below solving the equations for the angle a, s is the sigmoid_threshold:
% k0 and a0 are the optimized results for the sig slope and threshold
%
% s = 1/(1+exp(-k0(a-a0))
% 1+exp((-k0(a-a0)) = 1/s
% exp(-k0(a-a0)) = 1/s - 1
% -k0(a-a0) = log(1/s - 1)
% a-a0 = log(1/s - 1)/(-k0)
% a = log(1/s - 1)/(-k0) + a0
%
% note that in optimized sigmoid results the order of parameters is
% [amplitude slope threshold] (see sigmoid.m)
%    1         2       3       index

angle_rise = log( 1/sigmoid_threshold -1 ) / ( - opt_rise_params(2) ) + opt_rise_params(3);
angle_fall = log( 1/sigmoid_threshold -1 ) / ( - opt_fall_params(2) ) + opt_fall_params(3);

disp('The implied angles are those that segregate breath spiking activity (includes peak) and a quiet region (trough)')
disp(['For sigmoid_threshold = ' num2str(sigmoid_threshold) ' implied rise angle = ' num2str(angle_rise) ',  fall angle = ' num2str(angle_fall)]);
breath_angle_range = (angle_fall - angle_rise)/2;
breath_angle = (angle_fall + angle_rise)/2;
disp('For within the active spiking region of the breath cycle:')
disp(['breath angle = ' num2str(breath_angle) ', and  (plus/minus) breath_angle_range = ' num2str(breath_angle_range)]);

%20140923 if the peak occurs late in the breath angle like 5 radians then
% Shaina says these quantities are not calculated properly:
% assume that the trough does not overlap the 0 and 2pi borders
% find the left and right points of the trough
left_trough_angle = breath_angle+breath_angle_range;
right_trough_angle = breath_angle-breath_angle_range;
% however the right trough point needs to be moved from a negative angle to
% a positive one
right_trough_angle = mod(right_trough_angle, 2*pi);
trough_breath_angle = mean([left_trough_angle right_trough_angle]);
trough_breath_angle_range = (right_trough_angle - left_trough_angle)/2;
disp('For the trough region of the breath cycle:')
disp(['trough breath angle = ' num2str(trough_breath_angle) ', and  (plus/minus) trough_breath_angle_range = ' num2str(trough_breath_angle_range)]);

figure(multi_plot)
hold on
subplot(4,2,2)
angle_positions = mod([breath_angle-breath_angle_range breath_angle breath_angle+breath_angle_range], 2*pi);
plot(angle_positions, sigmoid_threshold*max(opt_rise_params(1),opt_fall_params(1)).*ones(1,3),'bo')
plot(angle_positions, sigmoid_threshold*max(opt_rise_params(1),opt_fall_params(1)).*ones(1,3),'p')

plot(trough_breath_angle, sigmoid_threshold*max(opt_rise_params(1),opt_fall_params(1)), 'co')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate plots of the breath angle ranges for the Fitness Image Video
frames=20
angle_list=[0:2*pi/frames:2*pi]+2*pi/(2*frames);
angle_range=1.5
for figures=1:1
    h=figure;
    hold on
    angle_positions = mod([breath_angle-breath_angle_range breath_angle breath_angle+breath_angle_range], 2*pi);
    plot(breath_angle_theta, smoothed_est_resp_firing_rate, 'linewidth',4)
    
    height=(max(smoothed_est_resp_firing_rate))*.65;
    first=angle_list(figures)-angle_range;
    second=angle_list(figures)+angle_range;
    
    if first<0
        plot([0 second],[height height],'k','linewidth',4)
        plot([mod(first,2*pi) 2*pi],[height height],'k','linewidth',4)
    else if second >2*pi
            plot([0 mod(second,2*pi)],[height height],'k','linewidth',4)
            plot([first 2*pi],[height height],'k','linewidth',4)
        else
            plot([first second],[height height],'k','linewidth',4)
        end
    end
    plot(angle_list(figures), height, 'r.','MarkerSize',30)
    title('Breath Activity', 'fontsize',16)
    set(gca, 'FontSize', 14)
    xlabel('Breath Angle (Radians)','fontsize',16)
    ylabel('Averaged Firing Rate (Hz)','fontsize',16)
    %ax=axis;
    %axis([0     7   0   110]);
    
    saveas(gcf, ['data/sort_breath/frames/angle_range_frames/angle_range_frames_' num2str(figures) '.png'], 'png')
    saveas(gcf, ['data/sort_breath/frames/angle_range_frames/angle_range_frames_' num2str(figures) '.eps'], 'eps')
end


fid=fopen('data/shared_file_pointers.dat','r');
line1=fgets(fid);
%cmd = ['load ' line1(1:end-1) ';'];
%eval(cmd);
%disp(['got tank from ' cmd]);
filename = line1(6:end -1);
%disp(['set filename to ' filename ]);
cellname = filename(1:end-4);

cmd=['fid=fopen(''../../static_data_summary.dat'',''a'');'];
eval(cmd)
phase_shift_cont=mod(peak_exp_angle-peak_control_angle+pi, 2*pi)-pi;
phase_shift=mod(peak_spont_angle-peak_control_angle+pi, 2*pi)-pi;
   fprintf(fid, '%s %g %g %g %g %g %g %g %g \n', cellname, peak_control_angle, peak_control_fr, peak_exp_angle, peak_exp_fr, phase_shift, peak_spont_angle, peak_spont_fr, phase_shift_cont);
   fclose(fid);
 

