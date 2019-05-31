% net_type_adder.m
% This script adds a line to a range X=startX, startX+1, ..., endX-1, endX of
% run_X/parameters.hoc files.  The line it adds identifies which type of
% NN the run had.  This program is an after the fact addition; the lines
% will subsequently be added as the networks are created by batch_run.py
% once I update batch_run.py to do this.  The added lines look like
% "// net_type pg_net" or "// net type gc_net" or if there is a special
% parameter choice of nc[15][:]=0 it could for example look like
% "// net_type pg_net0nc15". Don't forget to hand select folders before
% running this script.

% hand select below before starting program

absolute_path='/Users/morse/Documents/projects/VerhagenLab/theory/nrn/sim_results/20150710/';
sterotyped_folder='run_*';

startX=0;
endX=89;
startX=90;
endX=179;

%added_line=['// net_type ' 'pg_net0nc15'];
added_line=['// net_type ' 'gc_net'];
%added_line=['// net_type ' 'pg_net'];

% end hand selection section

folders=dir([absolute_path sterotyped_folder]);

for folder_index=1:length(folders)
    folder_name = folders(folder_index).name;
    % extract number from file name
    run_number=eval(strrep(folder_name,'run_','')); % folder_name looks like run_X
    if startX<=run_number && run_number<=endX        
        file_name = [absolute_path folder_name '/parameters.hoc'];
        fid=fopen(file_name,'a');
        fprintf(fid,[added_line '\n']);
        fclose(fid);
        % disp(['modified ' file_name]);
    end
end
