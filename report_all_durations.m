% report_all_durations.m
% read all the filenames of type tdt2mat_number.mat and report all the
% durations found in them one by one in a file
% all_durations.txt

% read all the filenames

%% list_of_files=dir('data/tdt2mat_data*.mat');

list_of_files=dir('data/20140714C1StaticR1.mat');
%% list_of_files=dir('data/201407*.mat');

% loop over them
clear intended_dur_cell; % cell matrix{file number} entries unique duration value, number of hits
% intended_dur_cellcell{filename_index} = [ [duration_value number_of_hits]; ...
% ] contents make it a rows x 2 columns matrix

for filename_index = 1:length(list_of_files);
    % get their durations
    [ actual_durations, unique_actual_durations, intended_durations] = duration_reporter( list_of_files(filename_index).name );
    rounded_act_durs = round(actual_durations);
    fid = fopen('all_durations.txt','at');
    fprintf(fid,[list_of_files(filename_index).name ': unique durations (within 1 ms), occurences\n']);
    for i=1: length(unique_actual_durations)
      ind_values = find(rounded_act_durs==unique_actual_durations(i));
      fprintf(fid, '%g, %g\n', unique_actual_durations(i), length(ind_values));
      intended_dur_cell{filename_index}(i,1) = unique_actual_durations(i);
      intended_dur_cell{filename_index}(i,2) =  length(ind_values);
    end
    fclose(fid);
    
end

figure;
for filename_index = 1:length(list_of_files)
    subplot(length(list_of_files),1,filename_index)
    hold on
    plot(intended_dur_cell{filename_index}(:,1),intended_dur_cell{filename_index}(:,2),'ro')
    plot(intended_dur_cell{filename_index}(:,1),intended_dur_cell{filename_index}(:,2),'b')
end

% now analyze the intended_dur_cell matricies to find suggested
% intended durations: 
clear intended_durs;  % look for the peaks
for filename_index = 1:length(list_of_files);
    % find the value of the local max
    local_max = max(intended_dur_cell{filename_index}(:,2));
    % everytime local_max is exceeded look for where drops back down again
    % then label that as a peak
    half_local_max = 0.5 * local_max;
    dur_deriv = diff(intended_dur_cell{filename_index}(:,2));
    % march through derivative variable dur_deriv: when positive on one
    % side and negative on another and if value is > half_local_max then
    % declare this as a possible peak
    file_specific_dur_index=1;
    for i=2: length(dur_deriv)
        if dur_deriv(i-1)>0 && dur_deriv(i)<0 && intended_dur_cell{filename_index}(i,2) > half_local_max
            intended_durs{filename_index}(1, file_specific_dur_index) = intended_dur_cell{filename_index}(i,1);
            file_specific_dur_index = file_specific_dur_index + 1;
        end
    end
    disp(['for ' list_of_files(filename_index).name ' the predicted intended durs are' char(10) '[' num2str(intended_durs{filename_index}) ']'])
end

