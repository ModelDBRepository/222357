% polar_plot2.m
% will read glomon and then look at the fourth and last column of the
% data/sortbreath/whatever.dat file to plot these values in a polar plot
% The radii are scaled arbitrarily to go from half the number of glomon's
% to the total number of glomons.  This arbitrary radii selection allows
% the stimuli's angles to be visible.  The fourth column represents the
% angle of S_ON + 50 ms and the last column represent the angle of S_ON +
% 150 + 50 ms.

load glomon.dat;
load data/sort_breath/20140616C1RandR3_20150609Epolar.dat;
a=X20140616C1RandR3_20150609Epolar ;
stimids = a(:,1);
start_angles=a(:,4);
stop_angles=a(:,end);


% just plot stimuli that start in this range:
%breath_angle=(13/24)*(2*pi)%(7/12)*(2*pi)
%breath_range=(2*pi)/6%(2*pi)/12
breath_angle=(1/12)*(2*pi); %+eps
breath_range=(2*pi)/12%(2*pi)/12

% polar plot of stimulus window angles
radii_index=1;
number_of_arcs=0;
max_r=1; % will be reset
min_r=1;
delta_r=1;
for run=1:2
    % first time through count the number of arcs that will be plotted
    % second time through (run=2): plot them
    if run==2
        arb_num=100 %length(glomon); % number_or_arcs;
        max_r=arb_num;
        min_r=max_r/3; % 0; max_r*.38;
        %delta_r = (max_r - min_r)/(number_of_arcs-1);
        delta_r = (max_r - min_r)/(arb_num-1);
        figure
        hold off
    end
    
    for index=1:length(glomon)
        angle_index=find(glomon(index)==stimids); % finds stimid in glomon in the big angle matrix
        if length(angle_index)
            start_angle=start_angles(angle_index);
            if test_breath_angle(start_angle, breath_angle, breath_range)
                stop_angle=stop_angles(angle_index);
                if start_angle>stop_angle
                    % goes through 2pi
                    stop_angle=stop_angle+2*pi;  % restores stop angle to original larger value before modulus
                end
                thetas=start_angle:(2*pi)/50:stop_angle;  % arbitrary resolution of 1/50th of a revolution
                r=max_r-(radii_index-1)*delta_r; % start at r=length(glomon) and move inwards
                radii_index = radii_index + 1;
                rs=r.*ones(1, length(thetas));
                
                if run==2
                    polar(thetas, rs, 'k')                
                    hold on
                else
                    number_or_arcs=number_of_arcs+1;
                end
                
            end
        else
            disp(['could not find stimid ' num2str(glomon(index)) ' in polar matrix'])
        end
        
    end
disp(['number of arcs ' num2str(number_of_arcs)]);    
end

axis equal
