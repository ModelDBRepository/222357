% batch_make_bi_polars.m
% sets up tanks for batch processing
% of breath increment (bi) polar plots

% absolute_path='/Users/morse/Documents/projects/VerhagenLab/20150618/louise_results/batch_runs/';
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/20150618/louise_results/20150623/batch_runs/';
% absolute_path='/Users/morse/Documents/projects/VerhagenLab/theory/nrn/sim_results/';
% date_path='20150708/';
% date_path='20150710/';
% date_path='20150714/';
% date_path='20150713a/';

% % absolute_path=['/Users/morse/Documents/projects/VerhagenLab/theory/nrn/sim_results/' date_path];
% absolute_path=['/home/tmm46/projects/VerhagenLab/' date_path '/batch_runs/'];

%date_path='20150716/'; % this is used in polar_plot2saver.m
%date_path='20150720/'; % this is used in polar_plot2saver.m
%date_path='control_incr5/'; % this is used in polar_plot2saver.m
%date_path='control_20_760_20/'; % this is used in polar_plot2saver.m
%date_path='cols_50_400'; % this is used in polar_plot2saver.m
%date_path='light2_added';
%%date_path='right1_added';
%date_path='soght1_added'; % for subset1-4 of light1_added
%date_path='cols_equal'; % for subset1-4 of light1_added
%date_path='cols_100_150_200_gc4';
%date_path='NSG_4xgc_light2_incr20';
%date_path='NSG_variable_xgc_cols';
%date_path='light2_incr20_4xgc_6';
%date_path='cols_4xgc_3';
%date_path='xgc';
%date_path='light2_reduced_subset1'; % light2 reduced inhib to m2
%date_path='light2_col_test';
%date_path = 'cols_s0';
%date_path = 'noinhib_cols';
%date_path='cols17to24';
%date_path='col1to6gc10to40light2_';
%date_path='gc70_90_light1';
%date_path='gc30_50_light1';
%date_path='eventrate2spikes';
%date_path='light2';
%date_path='tight2_incr20_added'; % tight2 instead of light2 to let other one finish unimpeded
%date_path='light1_incr5';
%%date_path='light1_recr5'; % to recreate increment 5 to double check for single_stim_elements etc collisions
%%date_path='light2_incr5_added';
%date_path='variable_xgc_cols_light1_redo';
%date_path='X5Hz7_5hw';
%date_path='X10Hz7_5hw_rw25';
%date_path='X10Hz30hw_rw25';
%date_path='X10Hz50hw';
%date_path='B100_200_300_500_s0_30gc_longer';
%date_path='d20160901';
%date_path='d20160909';
%date_path='d20160910';
%date_path='d20160929';
date_path='d20161008';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150716/full_net_50_400_100_300/batch_runs/';
%absolute_path='/home2/sms296/20150720/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/batch_runs/';
%absolute_path='/home2/sms296/20150720/control_incr5/batch_runs/';
%absolute_path='/home2/sms296/20150720/control_20_760_20/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_50_400/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_added/batch_runs/';
%absolute_path='/home2/sms296/20150720/light1_added/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_added/batch_runs/';
%%absolute_path='/home2/sms296/20150720/batch_runs/';
%absolute_path='/home2/sms296/20150810/light2_incr5/batch_runs/';
%%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light1_added/batch_runs/';
%%absolute_path='/home2/sms296/20150810/light2_incr5/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light1_added/batch_runs/subset3/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_equal/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_100_150_200/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/light2_4xgc/NSG_4xgc_light2/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2/batch_runs/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/NSG_4xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/light1_10_20_30_xgc/NSG_variable_xgc_cols/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150903/shaina_requests/light2_incr20_4xgc/light2_incr20_4xgc/subset6/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150903/shaina_requests/cols_4xgc/NSG_col_4xgc/subset3/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150903/variable_xgc_strength/NSG_variable_xgc_cols/'; % prev subsets
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150904eventrate2spikes/eventrate2spikes/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_reduced_inhib_to_m2/light2_reduced_inhib_to_m2/subset1/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/light2_col_test/NSG_col_4xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols_s0/NSG_col_S0_pg_gc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/noinhib_cols/NSG_col_4xgc_noinhib/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/cols17to24/NSG_col_17_through_23/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/gc10_to_40_cols1to6_light2/NSG_col1to6_gc10to40_light2/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/gc70_90_light1/light1_incr20_70and90xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/20150720/gc30_50_light1/light1_incr20_30and50xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20150902/light1_10_20_30_xgc/NSG_variable_xgc_cols/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/5Hz7.5HW/5Hz7.5hw/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/10Hz7.5HW/10Hz7.5hw/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151029/light1_cntrl_B100_to_500_s0_30xgc/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151031/light1_cntrl_B100_to_500_s0_30xgc_longer/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160826/NSG_ET_4xgc_light1/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160909/NSG_ET_4xgc_light1/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160910/NSG_ET_4xgc_light1/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/10Hz30HW/10Hz30hw/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20151015/10Hz50HW/10Hz50hw/';
%absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20160929/NSG_pg_no_lat_light1/';
absolute_path='/home/tmm46/projects/VerhagenLab/from_NSG/20161008/pg_long_no_lat/NSG_pg_no_lat_long/';

sterotyped_folder='run_*';

folders=dir([absolute_path sterotyped_folder]);
% I couldn't figure out how to get matlab to have a wildcard in the middle
% of the path to a file so will use the folder run_* to find how many files
% and just assume that all the tanks are there

% set the tank names in shared_file_pointers.dat

% 935 was B=0 S=0
big_file_initialized=0;
big_breath=[];
for folder_index=1:length(folders) % for some reason 933 choked
    if folders(folder_index).isdir
        folder_name=folders(folder_index).name;
        run_number = folder_name(5:end);
        disp(['Processing model in folder: ' folder_name ' corresponding to run number ' run_number]);
        [status, message, messageid]=mkdir(['data/' date_path]); % create subfolder if necessary
        fid=fopen(['data/' date_path '/shared_file_pointers.dat'],'w'); % date_path presence allows parallel use of this program
        tank_address = [absolute_path folder_name '/tdt2mat_data/0tdt2mat_data.mat'];
        fprintf(fid, [tank_address '\n\n']);
        fclose(fid)
        % tack on shared breath into tank
        if big_file_initialized
        else
            
%            load /home2/sms296/shared/_d_streams_d_BRTH_d_data.dat.orig;
             load /home/tmm46/shared/_d_streams_d_BRTH_d_data.dat
            big_file_initialized= 1;
%            big_breath=X_d_streams_d_BRTH_d_data_dat;
            big_breath=X_d_streams_d_BRTH_d_data;
        end
        cmd = ['load ' tank_address];
        eval(cmd)
        tdt2mat_data.streams.BRTH.data = big_breath; 
        
        % call analysis programs
        
        % these were the ones for tce polar plots
        %     evoke_act_half_2015E
        %     calc_peak_and_valley_fr_angles
        %Breath_angle_response_static_var_dur
        % load the type of network and B and S values
        load_type_of_net_B_S
        % use these for breath_increments method:
        if strcmp(B,'0') && strcmp(S,'0')
            disp('B=0 and S=0 condition detected so analysis skipped')
        else
            Assign_Ang_FR_Rand_June2015
            BASTATIC
            two_methods=2; % 1 alerts polar_plot1saver that calc_... was called
            % use two_methods=2 above for breath increments method
            polar_plot2saver
        end
        fclose('all')
        
    end
end
