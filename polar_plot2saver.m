% polar_plot2saver.m
% the revisions of this script over polar_plot1saver include the ability to
% save raw and derived data that goes into the polar plots into a folder
% hierarchy polar_plot_data/gc_net, polar_plot_data/pg_net
% See below for more documentation.
%
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

%polar([0],[max_radius],'w') % graphing the first point, even if white, determines the scale of the polar plot
% it was determined that it is convenient to have all the graphs at the
% same scale, say with a max radius of 40:
polar([0],[40],'w') % graphing the first point, even if white, determines the scale of the polar plot
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

% store data associated with graphs:
% Store all polar plot data in a new folder called polar_plot_data.
% Names of matrices:
% For polar plot raw data incorporate the indices into the names of the files:
% 
% bXsY.dat evoked radii
% bXsY_angle.dat at these angles
% bXsY_ste.dat evoked radii ste
% bXsY_ctrl.dat control radii at same above angles
% bXsY_ctrl_ste.dat control radii ste
% 
% For polar plot derived data also incorporate the indices into the names of the files:
% 
% bXsY_evoked_max_FR.dat
% bXsY_evoked_ave_r.dat
% more below

polar_data = ['polar_plot_data']; % date_path set in batch_make_pi_polars.m
filename = [polar_data '/' date_path '/' type_of_net '/' bXsY '.dat'];

% check existence of folders and if any don't exist then create
newSubFolder = [polar_data '/' date_path];
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
newSubFolder = [polar_data '/' date_path '/' type_of_net];
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
newSubFolder = ['images/' date_path];
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
save(strcat(polar_data, '/', date_path, '/',  type_of_net, '/', bXsY, '.dat'),'trace','-ascii');
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_angle.dat'),'trace_angle','-ascii');
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_ste.dat'),'trace_ste','-ascii');
evoked_max_FR = max(trace);
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_evoked_max_FR.dat'),'evoked_max_FR','-ascii');
evoked_ave_r = mean(trace);
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_evoked_ave_r.dat'),'evoked_ave_r','-ascii');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evoked

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate centroid
[centx, centy] = compute_centroid(trace_angle, trace);

plot([0 centx],[0 centy],'r','linewidth',4)

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

Control_directivity=directivity(trace_control);
Evoked_directivity=directivity(trace);

clear i; % reassigns i to square root of minus 1
Z=centx+centy*1i;
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

% bXsY_evoked_centroid_mag.dat
% bXsY_evoked_centroid_ang.dat
evoked_centroid_mag=abs(Z);
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_evoked_centroid_mag.dat'),'evoked_centroid_mag','-ascii');
evoked_centroid_ang=phi;
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_evoked_centroid_ang.dat'),'evoked_centroid_ang','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now for control

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z=centx_control+centy_control*1i;
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

% bXsY_ctrl_max_FR.dat
% bXsY_ctrl_ave_r.dat
ctrl_max_FR = max(trace_control);
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_ctrl_max_FR.dat'),'ctrl_max_FR','-ascii');
ctrl_ave_r = mean(trace);
save(strcat(polar_data, '/', date_path, '/', type_of_net, '/', bXsY, '_ctrl_ave_r.dat'),'ctrl_ave_r','-ascii');


% folder_index-1 is used below to keep figures syncronized with run_X
%savefig(['images/model_polar_plot_' run_number '_' num2str(two_methods) '.fig']); % save as fig
saveas(fig_handle,['images/' date_path '/' type_of_net '_' bXsY '_model_polar_plot_' run_number '_' num2str(two_methods) '.fig']); % save as fig
saveas(fig_handle,['images/' date_path '/' type_of_net '_' bXsY '_model_polar_plot_' run_number '_' num2str(two_methods)],'epsc'); % eps extension auto
% bXsY_ctrl_centroid_mag.dat
% bXsY_ctrl_centroid_ang.dat
ctrl_centroid_mag=abs(Z);
save([polar_data '/' date_path '/' type_of_net '/' bXsY '_ctrl_centroid_mag.dat'],'ctrl_centroid_mag','-ascii');
ctrl_centroid_ang=phi;
save([polar_data '/' date_path '/' type_of_net '/' bXsY '_ctrl_centroid_ang.dat'],'ctrl_centroid_ang','-ascii');

save([polar_data '/' date_path '/' type_of_net '/' bXsY '_ctrl.dat'],'trace_control','-ascii');
save([polar_data '/' date_path '/' type_of_net '/' bXsY '_ctrl_ste.dat'],'trace_control_ste','-ascii');
