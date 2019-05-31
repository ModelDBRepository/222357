% directivity_plot.m
% plots directivity 3D plot where the height (directivity) is a function of
% x and y (breath and light1 _peak_rates)

% a couple of choices below as to which date_NN_path is choosen
% date_NN_path = '20150708/gc_net/';
% %date_NN_path = '20150708/pg_net/'; % B, S > = < 50, 400
% data_title = 'B, S > = < 50, 400';
% 
% %date_NN_path = '20150709/gc_net/';
% date_NN_path = '20150709/pg_net0nc15/'; % B, S > = < 50, 400
% data_title = 'B, S > = < 50, 400';
% 
% date_NN_path = '20150710/gc_net/';
% date_NN_path = '20150710/pg_net/';
% data_title = 'B={60, 120, 180, 240, 300}, S={10, 20, ... , B-10, B}';
% 
% date_NN_path= '20150714/gc_net/';
% date_NN_path= '20150714/pg_net/';
% data_title = 'B, S = combinations of (50,400) and (100,300)'
%absolute_path=['/Users/morse/Documents/projects/VerhagenLab/theory/nrn/matlab_analysis/polar_plot_data/' date_NN_path];

date_NN_path= '20150713a/pg_net/'; % 20150713a corresponds to 20150714 on mac
absolute_path=['/home/tmm46/projects/VerhagenLab/theory/nrn/matlab_analysis/polar_plot_data/' date_NN_path];

files = dir(absolute_path);

for index=3:length(files)  % first two are . and ..
    cmd=['load ' absolute_path files(index).name ';'];
    eval(cmd)
end

% the B, S values (or there transformation to number of average events)
% will be the x,y coordinates

% find all the B, S values present in the data

B_values_str={}; % first create string cells of the B and S values in the file names
S_values_str={};
B_values=[];
S_values=[];
for index=3:length(files)
    file_name=files(index).name;
    % from http://www.mathworks.com/matlabcentral/newsreader/view_thread/298630
    % B values:
    expr1='b';
    expr2='s';
    expr = [expr1 '(.*?)' expr2];
    tok = regexp(file_name,expr,'tokens');
    B_values = [B_values; eval(char(tok{:}))];
    % S values:
    expr1='s';
    expr2='_';
    expr = [expr1 '(.*?)' expr2];
    tok = regexp([file_name '_'],expr,'tokens'); % _ added to end in case none
    S_values = [S_values; eval(strrep(char(tok{1}),'.dat',''))]; % when
    % there is no _  there is actually a .dat at the end of the string that
    % needs to be removed
end
B_values=unique(B_values)';
S_values=unique(S_values)';

% find the directivity for each pair of these B, S's
B_X=[]; % the x values for the scatter plot
S_Y=[]; % the y values for the scatter plot
D_Z=[]; % the z values for the scatter plot
D_Z_ctrl=[]; % the z values for the scatter plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evoked and Control directivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % below is special case for B=60,120,180,240,300 and S=10,20,...,B-10, B
% for each B:
% for B_index=1:length(B_values)
%     B_tmp = B_values(B_index);
%     %    for S_index=1:length(S_values)
%     for S_index=1:B_tmp/10 % makes S_tmp = {10, 20, ..., B_tmp-10, B_tmp}
%         S_tmp = S_values(S_index);
%         cmd = ['tmp_max_FR = b' num2str(B_tmp) 's' num2str(S_tmp) '_evoked_max_FR;'];
%         eval(cmd);
%         cmd = ['tmp_ave_FR = b' num2str(B_tmp) 's' num2str(S_tmp) '_evoked_ave_r;'];
%         eval(cmd);
%         directivity_tmp = tmp_max_FR*tmp_max_FR/tmp_ave_FR;
%
%         cmd = ['tmp_max_FR_ctrl = b' num2str(B_tmp) 's' num2str(S_tmp) '_ctrl_max_FR;'];
%         eval(cmd);
%         cmd = ['tmp_ave_FR_ctrl = b' num2str(B_tmp) 's' num2str(S_tmp) '_ctrl_ave_r;'];
%         eval(cmd);
%         directivity_tmp_ctrl = tmp_max_FR_ctrl*tmp_max_FR_ctrl/tmp_ave_FR_ctrl;
%
%         B_X = [B_X; B_tmp];
%         S_Y = [S_Y; S_tmp];
%         D_Z = [D_Z; directivity_tmp];
%         D_Z_ctrl = [D_Z_ctrl; directivity_tmp_ctrl];
%     end
% end

% below for special case of B,S = (50, 50), (50, 400), (400, 50),
% (400,400), (100, 100), (100, 300), (300, 100), (300, 300)
%
B_values = [50 50 400 400 100 100 300 300];
S_values = [50 400 50 400 100 300 100 300];

for both_index=1:length(B_values) % index for both S and B
    B_tmp = B_values(both_index);
    %    for S_index=1:length(S_values)
    S_tmp = S_values(both_index);
    cmd = ['tmp_max_FR = b' num2str(B_tmp) 's' num2str(S_tmp) '_evoked_max_FR;'];
    eval(cmd);
    cmd = ['tmp_ave_FR = b' num2str(B_tmp) 's' num2str(S_tmp) '_evoked_ave_r;'];
    eval(cmd);
    directivity_tmp = tmp_max_FR*tmp_max_FR/tmp_ave_FR;
    
    cmd = ['tmp_max_FR_ctrl = b' num2str(B_tmp) 's' num2str(S_tmp) '_ctrl_max_FR;'];
    eval(cmd);
    cmd = ['tmp_ave_FR_ctrl = b' num2str(B_tmp) 's' num2str(S_tmp) '_ctrl_ave_r;'];
    eval(cmd);
    directivity_tmp_ctrl = tmp_max_FR_ctrl*tmp_max_FR_ctrl/tmp_ave_FR_ctrl;
    
    B_X = [B_X; B_tmp];
    S_Y = [S_Y; S_tmp];
    D_Z = [D_Z; directivity_tmp];
    D_Z_ctrl = [D_Z_ctrl; directivity_tmp_ctrl];
end

% Each graph separately and together in one figure
% evoked
fig_handle=figure;
subplot(3,1,1)
scatter3(B_X, S_Y, D_Z)
hold on
xlabel('breath\_peak\_rate')
ylabel('light1\_peak\_rate')
zlabel('Directivity')
title(['Evoked directivity, ' data_title ', ' strrep(date_NN_path,'_','\_')])

% control
subplot(3,1,2)
scatter3(B_X, S_Y, D_Z_ctrl)
hold on
xlabel('breath\_peak\_rate')
ylabel('light1\_peak\_rate')
zlabel('Directivity')
title(['Control directivity, ' data_title])

% both
subplot(3,1,3)
s = repmat([20],numel(B_X),1)';
c = repmat([1 0 0],numel(B_X),1);
hold off
scatter3(B_X, S_Y, D_Z, s, c)
hold on
%s = repmat([8],numel(B_X),1)';
% s = repmat([20],numel(B_X),1)';
c = repmat([0 0 0],numel(B_X),1);
scatter3(B_X, S_Y, D_Z_ctrl, s, c);
xlabel('breath\_peak\_rate')
ylabel('light1\_peak\_rate')
zlabel('Directivity')
title(['Evoked (red), Control (black) ' data_title]);

saveas(fig_handle,['images/directivity_plot_' strrep(date_NN_path,'/','_') '_evoked_and_ctrl.fig'],'fig'); % save as fig
saveas(fig_handle,['images/directivity_plot_' strrep(date_NN_path,'/','_') '_evoked_and_ctrl.eps'],'epsc'); % eps extension auto

% a new difference plot
fig_handle=figure;
s = repmat([20],numel(B_X),1)';
c = repmat([1 0 0],numel(B_X),1);
hold off
delta_D_Z = D_Z-D_Z_ctrl;
scatter3(B_X, S_Y, delta_D_Z, s, c)
hold on
xlabel('breath\_peak\_rate')
ylabel('light1\_peak\_rate')
zlabel('Delta Directivity')
title(['Evoked minus control directivity, ' data_title ', ' strrep(date_NN_path,'_','\_')])

saveas(fig_handle,['images/delta_directivity_plot_' strrep(date_NN_path,'/','_') '_evoked_minus_ctrl.fig'],'fig'); % save as fig
saveas(fig_handle,['images/delta_directivity_plot_' strrep(date_NN_path,'/','_') '_evoked_minus_ctrl.eps'],'epsc'); % eps extension auto

% difference plot as surface plot
fig_handle=figure;
[xx, yy]=meshgrid(B_X, S_Y);
[row_xx, col_xx]=size(xx);
[row_yy, col_yy]=size(yy);
zz=NaN.*zeros(size(xx)); % could have also had yy in place of xx
for row_index=1: row_xx
    for col_index=1:col_xx
        B_tmp=xx(row_index, col_index);
        S_tmp=yy(row_index, col_index);
        % is there a Z value for this B_tmp, S_tmp pair?
        for search_index=1: length(B_X)
            if (B_X(search_index)==B_tmp) && (S_Y(search_index)==S_tmp)
                zz(row_index, col_index)=delta_D_Z(search_index);
            end
        end
    end
end
hold off
surfl(xx, yy, zz)
hold on
colormap(pink)    % change color map
shading interp
xlabel('breath\_peak\_rate')
ylabel('light1\_peak\_rate')
zlabel('Delta Directivity')
title(['Evoked minus control directivity, ' data_title ', ' strrep(date_NN_path,'_','\_')])

saveas(fig_handle,['images/delta_directivity_surface_' strrep(date_NN_path,'/','_') '_evoked_minus_ctrl.fig'],'fig'); % save as fig
saveas(fig_handle,['images/delta_directivity_surface_' strrep(date_NN_path,'/','_') '_evoked_minus_ctrl.eps'],'epsc'); % eps extension auto

% 2d line plots
fig_handle=figure;
hold on
last_index=0;
% specific for the runs B=60, 120, 180, 240, 300; S=10, 20, ..., B-10, B
% 
% for B_index=1:length(B_values)
%     if B_index==1
%         plot(S_Y([1:B_values(B_index)/10]),D_Z(last_index+[1:B_values(B_index)/10]))
%     else
%         plot(S_Y(last_index+[1:B_values(B_index)/10]),D_Z(last_index+[1:B_values(B_index)/10]))
%     end
%     %disp('another set of indicies')
%     %[1:B_values(B_index)/10]
%     last_index = last_index + B_values(B_index)/10;
% end

% specific for the runs 
%B_values = [50 50 400 400 100 100 300 300];
%S_values = [50 400 50 400 100 300 100 300];
for both_index=2:2:length(B_values)
    plot(S_Y((both_index-1):both_index),D_Z((both_index-1):both_index));
end
% for both_index=1:length(B_values) % index for both S and B

xlabel('light1\_peak\_rate')
ylabel('Directivity')
title('Evoked Directivity - the breath is given by the value of the last light stimulus in the line')
saveas(fig_handle,['images/directivity_2D_' strrep(date_NN_path,'/','_') '_evoked.fig'],'fig'); % save as fig
saveas(fig_handle,['images/directivity_2D_' strrep(date_NN_path,'/','_') '_evoked.eps'],'epsc'); % eps extension auto

% 
% fig_handle=figure;
% hold on
% last_index=0;
% for B_index=1:length(B_values)
%     if B_index==1
%         plot(S_Y([1:B_values(B_index)/10]),D_Z_ctrl(last_index+[1:B_values(B_index)/10]))
%     else
%         plot(S_Y(last_index+[1:B_values(B_index)/10]),D_Z_ctrl(last_index+[1:B_values(B_index)/10]))
%     end
%     %disp('another set of indicies')
%     %[1:B_values(B_index)/10]
%     last_index = last_index + B_values(B_index)/10;
% end
% 
% xlabel('light1\_peak\_rate')
% ylabel('Directivity')
% title('Control Directivity - the breath is given by the value of the last light stimulus in the line')
% saveas(fig_handle,['images/directivity_2D_' strrep(date_NN_path,'/','_') '_ctrl.fig'],'fig'); % save as fig
% saveas(fig_handle,['images/directivity_2D_' strrep(date_NN_path,'/','_') '_ctrl.eps'],'epsc'); % eps extension auto
% 
% % delta
% fig_handle=figure;
% hold on
% last_index=0;
% for B_index=1:length(B_values)
%     if B_index==1
%         plot(S_Y([1:B_values(B_index)/10]),D_Z(last_index+[1:B_values(B_index)/10])-D_Z_ctrl(last_index+[1:B_values(B_index)/10]))
%     else
%         plot(S_Y(last_index+[1:B_values(B_index)/10]),D_Z(last_index+[1:B_values(B_index)/10])-D_Z_ctrl(last_index+[1:B_values(B_index)/10]))
%     end
%     %disp('another set of indicies')
%     %[1:B_values(B_index)/10]
%     last_index = last_index + B_values(B_index)/10;
% end
% 
% xlabel('light1\_peak\_rate')
% ylabel('Delta directivity')
% title('Evoked minus Control Directivity - the breath is given by the value of the last light stimulus in the line')
% saveas(fig_handle,['images/delta_directivity_2D_' strrep(date_NN_path,'/','_') '_evoked_minus_ctrl.fig'],'fig'); % save as fig
% saveas(fig_handle,['images/delta_directivity_2D_' strrep(date_NN_path,'/','_') '_evoked_minus_ctrl.eps'],'epsc'); % eps extension auto
% 
% 
