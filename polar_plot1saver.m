% polar_plot1.m
% experimental breath sorted instantaneous firing rate
% and control breath sorted instantaneous firing rate polar plots
fig_handle=figure;

% use method to match program to desired polar graph (resp, contr),(breath angle increments, breath_angle)
% two_methods = 1; % change this method for below choices
% *** set 'method' variable prior to calling this script ***

switch two_methods
    case 1
    % if using calc_peak... use these:
    trace_angle = breath_angle_theta;
    trace_control = smoothed_est_contr_firing_rate'; % or smoothed_est_resp_firing_rate
    %
    % if using calc_peak... use these:
    % trace_angle = breath_angle_theta;
    trace = smoothed_est_resp_firing_rate'; % or smoothed_est_resp_firing_rate
    case 2
    % if using Breath_increments use these
    trace_angle = xcenters;
    trace_control = FR_breath_angles_control; % or trace_angle_control
    trace_control_ste=ste_error_control;
    %
    % if using Breath_increments use these
    % trace_angle = xcenters;
    trace = FR_breath_angles; % or FR_breath_angles_control
    trace_ste = ste_error_exp;
end

% first find the max radius
all_radii=[]; % find all big radii so a max can be plotted to scale graph
all_radii=[all_radii trace_control];
all_radii=[all_radii trace];
all_radii=[all_radii trace+trace_ste];
all_radii=[all_radii trace_control+trace_control_ste];

max_radius=max(all_radii);

polar([0],[max_radius],'w') % graphing the first point, even if white, determines the scale of the polar plot
hold on
%title ('Experimental Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
% polar([trace trace(1)], [trace_angle+ste_error_exp trace_angle(1)+ste_error_exp(1)],'y')

polar([trace_angle trace_angle(1)], [trace trace(1)],'r')

h=polar([trace_angle trace_angle(1)], [trace+trace_ste trace(1)+trace_ste(1)]);
h.Color=[1 0 0];
h=polar([trace_angle trace_angle(1)], [trace-trace_ste trace(1)-trace_ste(1)]);
h.Color=[1 0 0];

hold on
% polar([trace trace(1)], [trace_angle-ste_error_exp trace_angle(1)-ste_error_exp(1)],'y')
axis equal
% calculate centroid
[centx, centy] = compute_centroid(trace_angle, trace);

plot([0 centx],[0 centy],'r','linewidth',4);

% Control

%title('Control Breath Sorted Instantaneous Firing Rate','fontsize',label_fontsize)
% polar([trace trace(1)], [trace_angle_control+ste_error_control trace_angle_control(1)+ste_error_control(1)], 'y')

hold on
% trace(1) repeats the first value to close the polar plot
polar([trace_angle trace_angle(1)],[trace_control trace_control(1)],'k')
h=polar([trace_angle trace_angle(1)], [trace_control+trace_control_ste trace_control(1)+trace_control_ste(1)]);
nblack = (.5/3)^(0.5); % near black: fraction from black magnitude from nblack being in each color coordinate
h.Color=[0 0 0]; %[nblack nblack nblack];
h=polar([trace_angle trace_angle(1)], [trace_control-trace_control_ste trace_control(1)-trace_control_ste(1)]);
h.Color=[0 0 0] ; % nblack nblack nblack];

% polar([trace trace(1)], [trace_angle_control-ste_error_control trace_angle_control(1)-ste_error_control(1)], 'y')
% ylabel('spikes/(num_of_stims*control_response_period) (Hz)', 'Interpreter', 'none','fontsize',label_fontsize)
% xlabel(['breath angle bins (radians) 2*pi/' num2str(num_of_breath_increments) ' bin widths','fontsize',label_fontsize])
[centx_control, centy_control]=compute_centroid(trace_angle, trace_control);
axis equal
plot([0 centx_control],[0 centy_control],'k','linewidth',4)
% compute_centroid(trace_angle, trace_control)

Control_directivity=directivity(trace_control)
Evoked_directivity=directivity(trace)

clear i; % reassigns i to square root of minus 1
Z=centx+centy*i;
clear angle; % makes the angle function available
phi=angle(Z); % returns angle of Z in -pi to pi
if phi<0
    phi = phi+2*pi;
    % returns pi to 0 <= phi < 2*phi
end
disp(['For response: ********************************'])
disp(['The angle of the centroid is ' num2str(phi)])
disp(['magnitude of the centroid is ' num2str(abs(Z))])
disp(['average polar radius is ' num2str(mean(trace))])
disp(['magnitude of the centroid/average polar radii is ' num2str(abs(Z)/mean(trace))])

% now for control
Z=centx_control+centy_control*i;
phi=angle(Z); % returns angle of Z in -pi to pi
if phi<0
    phi = phi+2*pi;
    % returns pi to 0 <= phi < 2*phi
end
disp(['For control: ********************************'])
disp(['The angle of the centroid is ' num2str(phi)])
disp(['magnitude of the centroid is ' num2str(abs(Z))])
disp(['average polar radius is ' num2str(mean(trace_control))])
disp(['magnitude of the centroid/average polar radii is ' num2str(abs(Z)/mean(trace_control))])

% folder_index-1 is used below to keep figures syncronized with run_X
savefig(['images/model_polar_plot_' run_number '_' num2str(two_methods) '.fig']); % save as fig
saveas(fig_handle,['images/model_polar_plot_' run_number '_' num2str(two_methods)],'epsc'); % eps extension auto

