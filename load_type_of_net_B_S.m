% load_type_of_net_B_S.m
% load the type of network (type_of_net) and B and S values
% all of these variables are strings

% because the way these are determined may be very specific to a batch
% simulation run they methods used to determine them may be very heuristic
filetext=fileread([absolute_path folder_name '/num_of_columns.hoc']);

% new: read in the number of columns and use that to form the type of net 
matlab_command=regexprep(filetext, '//',';%') % by substituting the comment
% character the hoc command is a perfectly good matlab command
eval(matlab_command); % now n holds the num_of_columns (extra columns actually)
% the total number of columns is this n + 1
total_num_of_columns = num2str(n+1);


% get the relevant parameter file information
filetext=fileread([absolute_path folder_name '/parameters.hoc']);

% now use search expressions to find things like the peak rates
expr = '[^\n]*peak_rate[^\n]*';

fileread_info = regexp(filetext,expr,'match');
eval([fileread_info{1} ';'])
eval([fileread_info{2} ';'])
eval([fileread_info{3} ';'])

% optional addition of searching for gc synapse strength (weight) sets
% gc_on if exists
expr = '[^\n]*gc_on[^\n]*';
gc_on_exists=0;
fileread_info = regexp(filetext,expr,'match');
% uncomment the below if desired to have gc g in file names
 if length(fileread_info)
     eval(strrep(fileread_info{1},'//','%'));
     gc_on_exists=1;
 end

% 
% now breath_peak_rate, light1_peak_rate, light2_peak_rate are available
B=num2str(breath_peak_rate);
% ************ IMPORTANT *********** slect the proper one of the below for
% light1 for light1, light2 for light2
S=num2str(light1_peak_rate);
%S=num2str(light2_peak_rate);
xgc='';
if gc_on_exists
    xgc=num2str(gc_on)
end
    % now for the network type

% use new "// net type " keyword
expr = '[^\n]*// net_type [^\n]*';
fileread_info = regexp(filetext,expr,'match');
% makes for example fileread_info = '// net_type gc_net'
type_of_net=char(strrep(fileread_info,'// net_type ',''));
% now assignment is for example 'gc_net'

% old method which did not work but may someday
if 0
    type_of_net = '_net'; % start with _net and then add descriptions
    
    expr = '[^\n]*toggle_gc_connection[^\n]*';
    
    fileread_info = regexp(filetext,expr,'match');
    
    if length(fileread_info)
        % if true then a pg network because gc toggled off
        type_of_net = ['pg' type_of_net];
        expr = '[^\n]*nc\[15\][^\n]*';
        fileread_info = regexp(filetext,expr,'match');
        if length(fileread_info)
            type_of_net = [type_of_net '0nc15'];
        end
    else
        % check to see if a gc network
        expr = '[^\n]*toggle_pg_connection[^\n]*';
        fileread_info = regexp(filetext,expr,'match');
        if length(fileread_info)
            type_of_net = ['gc' type_of_net];
        end
    end
end % if 0 commenting out of old method

% in general don't need column in net folder name because just 2 columns
% if do have columns varying comment in the below, otherwise leave commented:
type_of_net = [type_of_net '_' total_num_of_columns '_cols_gc_' xgc] % added to folder name

if gc_on_exists && 0 % 20150922 decided to put the gc in the folder name above
    bXsY=['b' B 's' S 'g' xgc]; % bXsYgZ is used to form filenames
else
    bXsY=['b' B 's' S]; % bXsY is used to form filenames
end
