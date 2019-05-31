% batch_make_tanks.m
% will put a tank in each folder when a stereotyped folder is passed to
% this function
% Usage:
% make_batch_tanks(absolute_path, stereotyped_folder)

%function batch_tanks(absolute_path, sterotyped_folder)
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/20150618/louise_results/batch_runs/';
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/20150618/louise_results/20150623/batch_runs/';
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/theory/nrn/sim_results/';
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/theory/nrn/sim_results/20150708/';
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/theory/nrn/sim_results/20150710/';
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/theory/nrn/sim_results/20150714/';
% absolute_path='/home/tmm46/projects/VerhagenLab/20150713a/batch_runs/';
% absolute_path='/home/tmm46/projects/VerhagenLab/20150716/full_net_50_400_100_300/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_added/batch_runs/';
% absolute_path='/home/tmm46/projects/VerhagenLab/20150720/batch_runs/';
% absolute_path='/home2/sms296/20150720/control_incr5/batch_runs/';
% absolute_path='/home2/sms296/20150720/control_20_760_20/batch_runs/';
% absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_50_400/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_300_100/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_400_50/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_100_300/batch_runs/';
%absolute_path='/home2/sms296/20150810/light1_incr5/batch_runs/subset1/';
%absolute_path='/home2/sms296/20150810/light2_incr5/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light1_added/batch_runs/subset4/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_equal/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_100_150_200/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/NSG_4xgc/';
% absolute_path='/home2/sms296/20150720/light1_added/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/light2_4xgc/NSG_4xgc_light2/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/light1_10_20_30_xgc/NSG_variable_xgc_cols/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150903/shaina_requests/light1_incr20_4xgc/light1_incr20_4xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150903/shaina_requests/light2_incr20_4xgc/light2_incr20_4xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150903/shaina_requests/cols_4xgc/NSG_col_4xgc/';
%absolute_path='/home/tmm46//projects/VerhagenLab/from_NSG/20150903/variable_xgc_strength/NSG_variable_xgc_cols/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150904eventrate2spikes/eventrate2spikes/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_cols/NSG_4xgc_light2_cols';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_reduced_inhib_to_m2/light2_reduced_inhib_to_m2/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_col_test/NSG_col_4xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_s0/NSG_col_S0_pg_gc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/noinhib_cols/NSG_col_4xgc_noinhib/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols17to24/NSG_col_17_through_23/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/gc10_to_40_cols1to6_light2/NSG_col1to6_gc10to40_light2/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/gc30_50_light1/light1_incr20_30and50xgc/';
% absolute_path='/home/tmm46/projects/VerhagenLab/20150720/gc70_90_light1/light1_incr20_70and90xgc/';
% absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/5Hz7.5HW/5Hz7.5hw/';
% absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/10Hz7.5HW/10Hz7.5hw/';
% absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/10Hz30HW/10Hz30hw/';
% absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/10Hz50HW/10Hz50hw/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151029light1_cntrl_B100_to_500_s0_30xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151031/light1_cntrl_B100_to_500_s0_30xgc_longer/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160826/NSG_ET_4xgc_light1/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160909/NSG_ET_4xgc_light1/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160910/NSG_ET_4xgc_light1/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160929/NSG_pg_no_lat_light1/';
absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20161008/pg_long_no_lat/NSG_pg_no_lat_long/';

sterotyped_folder='run_*';

folders=dir([absolute_path sterotyped_folder]);

for folder_index=1:length(folders)
    if folders(folder_index).isdir
        folder_to_process = [absolute_path folders(folder_index).name '/tdt2mat_data'];
        disp([' *** ' folder_to_process]);
        tdt2mat_data = nrn2tdt(folder_to_process,'no auto write');
        save([folder_to_process '/0tdt2mat_data.mat'], 'tdt2mat_data'); % put
        % number in tank file name to legitimize X in filename string in analysis programs
    end
end

