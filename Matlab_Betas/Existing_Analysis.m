function Existing_Analysis(hObject)
index_selected = get(hObject,'value'); 
Analysis = get(hObject,'string'); Analysis = char(Analysis(index_selected));
if strcmpi(Analysis,'No Analyses') == 1;
    return
end
path = ['C:\KnobAnalysis\ConfigFiles\' Analysis];
load(path);
devices = unique({config.data.device});
FinalDates = {config.FinalDates}; FinalDates = FinalDates{:};
All_Animals_Value = config.all_animals_value;
pos = get(0,'Screensize');                                              %Grab the screensize.
h = 5;                                                                 %Set the height of the figure, in centimeters.
w = 15;                                                                 %Set the width of the figure, in centimeters.
for d = 1:length(devices)
    fig = figure('numbertitle','off','units','centimeters',...
        'name',['Dexterity: Subject Overview'],'menubar','none',...
        'position',[pos(3)/2-w/2, pos(4)/2-h/2, w, h]);
    tgroup = uitabgroup('Parent', fig);
    tabs = {config.data.rat};
    plotdata = config.data;
    
    for i = 1:length(tabs);
        tab(i) = uitab('Parent', tgroup, 'Title', sprintf('Animal %s', tabs{i}));
        tabs_info(i,:) = get(tab(i));
        laststages(i) = plotdata(i).stage(length(plotdata(i).stage));
        numberofsessions(i) = length(plotdata(i).stage);
        AnimalName(i) = uicontrol('Parent', tab(i), 'Style', 'text', 'String', sprintf('Animal Name: %s', tabs{i}), ...
            'HorizontalAlignment', 'left', 'units', 'normalized', 'Position', [.05 .85 .4 .1],...
            'Fontsize', 10, 'Fontweight', 'normal') ;
        LastSessionRun(i) = uicontrol('Parent', tab(i), 'Style', 'text', 'String', sprintf('Last Session Run: %s', laststages{i}), ...
            'HorizontalAlignment', 'left', 'units', 'normalized', 'Position', [.05 .65 .6 .1],...
            'Fontsize', 10, 'Fontweight', 'normal') ;
        LastDateRun(i) = uicontrol('Parent', tab(i), 'Style', 'text', 'String', sprintf('Last Date Run: %s', FinalDates{i}), ...
            'HorizontalAlignment', 'left', 'units', 'normalized', 'Position', [.05 .45 .4 .1],...
            'Fontsize', 10, 'Fontweight', 'normal') ;
        NumberOfSessions(i) = uicontrol('Parent', tab(i), 'Style', 'text', 'String', sprintf('Number of Sessions: %d', numberofsessions(i)), ...
            'HorizontalAlignment', 'left', 'units', 'normalized', 'Position', [.05 .25 .4 .1],...
            'Fontsize', 10, 'Fontweight', 'normal') ;
        DevicesRuns(d) = uicontrol('Parent', tab(i), 'Style', 'text', 'String', sprintf('Devices: %s', devices{d}), ...
            'HorizontalAlignment', 'left', 'units', 'normalized', 'Position', [.05 .05 .4 .1],...
            'Fontsize', 10, 'Fontweight', 'normal') ;
        General_Analysis(i) = uicontrol('Parent', tab(i),'style','pushbutton','string','General Analysis','HorizontalAlignment', 'left',...
            'units','normalized','position',[.7 .625 .25 .25],'fontsize',10, 'callback', {@GeneralAnalysis,devices(d)},'userdata',plotdata(i));
    end
   if All_Animals_Value == 1;
        if length({config.data.rat}) == 1;
            uiwait(msgbox('Only one animal is loaded so a separate figure with all animals will not be displayed',...
                'Only One Animal Found'));
        else
            
            data = plotdata;
            pos = get(0,'Screensize');                                              %Grab the screensize.
            h = 10;                                                                 %Set the height of the figure, in centimeters.
            w = 15;                                                                 %Set the width of the figure, in centimeters.
%             for d = 1:length(devices)                                                   %Step through the devices.
                fig = figure('numbertitle','off','units','centimeters',...
                    'name','Dexterity: View All Subjects','menubar','none',...
                    'position',[pos(3)/2-w/2, pos(4)/2-h/2, w, h]);                     %Create a figure.
                ui_h = 0.07*h;                                                          %Set the heigh of the uicontrols.
                fontsize = 0.6*28.34*ui_h;                                              %Set the fontsize for all uicontrols.
                sp1 = 0.02*h;                                                           %Set the vertical spacing between axes and uicontrols.
                sp2 = 0.01*w;                                                           %Set the horizontal spacing between axes and uicontrols.
%                 pos = [7*sp2,3*sp1,w-8*sp2,h-ui_h-5*sp1];                               %Set the position of the axes.
                pos = [7*sp2,7*sp1,w-8*sp2,h-ui_h-10*sp1]; 
                ax = axes('units','centimeters','position',pos,'box','on',...
                    'linewidth',2);                                                     %Create axes for showing the log events histogram.
                obj = zeros(1,7);                                                       %Create a matrix to hold timescale uicontrol handles.
                str = {'Overall Hit Rate',...
                    'Total Trial Count',...
                    'Mean Peak Force',...
                    'Median Peak Force',...
                    'Trial Count',...
                    'Hits in First 5 Minutes',...
                    'Trials in First 5 Minutes',...
                    'Max. Hits in Any 5 Minutes',...
                    'Max. Trials in Any 5 Minutes',...
                    'Max. Hit Rate in Any 5 Minutes',...
                    'Min. Inter-Trial Interval (Smoothed)',...
                    'Mean Peak Impulse',...
                    'Median Peak Impulse',...
                    'Peak Velocity',...
                    'Latency to Hit'};                                             %List the available plots for the pull data.
                if any(strcmpi(devices{d},{'knob','lever'}))                            %If we're plotting knob data...
                    str(2:3) = {'Mean Peak Angle','Median Peak Angle'};                 %Set the plots to show "angle" instead of "force".
                elseif ~any(strcmpi(devices{d},{'knob','lever','pull'}))                %Otherwise, if this isn't pull, knob, or lever data...
                    str(2:3) = {'Mean Signal Peak','Median Signal Peak'};               %Set the plots to show "signal" instead of "peak force".
                end
                pos = [sp2, h-sp1-ui_h, 2*(w-6*sp2)/6, ui_h];                           %Set the position for the pop-up menu.
                obj(1) = uicontrol(fig,'style','popup','string',str,...
                    'units','centimeters','position',pos,'fontsize',fontsize);      %Create pushbuttons for selecting the timescale.
                str = {'Session','Day','Week','Export'};                            %List the timescale labels.
                for i = 2:5                                                             %Step through the 3 timescales.
                    pos = [i*sp2+i*(w-6*sp2)/6, h-sp1-ui_h, (w-6*sp2)/6, ui_h];         %Set the position for each pushbutton.
                    obj(i) = uicontrol(fig,'style','pushbutton','string',str{i-1},...
                        'units','centimeters','position',pos,'fontsize',fontsize);      %Create pushbuttons for selecting the timescale.
                end
                pos = [sp2, sp1, 2*(w-6*sp2)/6, ui_h]; 
                obj(6) = uicontrol(fig,'style','radiobutton','string','Concatenate Dates',...
                    'units','centimeters','position',pos,'fontsize',fontsize);
                pos = [30*sp2, sp1, 2*(w-6*sp2)/6, ui_h];
                obj(7) = uicontrol(fig,'style','radiobutton','string','Group Animals',...
                    'units','centimeters','position',pos,'fontsize',fontsize);
                set(obj(1),'callback',{@Set_Plot_Type,obj,data});                            %Set the callback for the pop-up menu.
                set(obj(2:4),'callback',{@Plot_Timeline,obj,[],data});                       %Set the callback for the timescale buttons.
                set(obj(5),'callback',{@Export_Data,ax,obj,data});                           %Set the callback for the export button.
                set(obj(6),'callback',{@TrainingDays,obj,data});
                set(obj(7),'callback',{@GroupAnimals,obj,data});
                set(fig,'userdata',data);                                                              %Save the plot data to the figure's 'UserData' property.
                Plot_Timeline(obj(2),[],obj,[],data);                                        %Call the function to plot the session data in the figure.
                set(fig,'ResizeFcn',{@Resize,ax,obj});
%             end
        end
    end
end

function TrainingDays(~,~,obj,data)
i = strcmpi(get(obj,'fontweight'),'bold');
Plot_Timeline(obj(i),[],obj,[],data);     

function GroupAnimals(~,~,obj,data)
i = strcmpi(get(obj,'fontweight'),'bold');
Plot_Timeline(obj(i),[],obj,[],data);

%% This function is called when the user selects the General Analysis Pushbutton
function GeneralAnalysis(hObject,~,devices)
data = hObject.UserData;
pos = get(0,'Screensize');                                              %Grab the screensize.
h = 10;                                                                 %Set the height of the figure, in centimeters.
w = 15;                                                                 %Set the width of the figure, in centimeters.
for d = 1:length(devices)                                                   %Step through the devices.
    fig = figure('numbertitle','off','units','centimeters',...
        'name','Dexterity: Single Subject','menubar','none',...
        'position',[pos(3)/2-w/2, pos(4)/2-h/2, w, h]);                     %Create a figure.
    ui_h = 0.07*h;                                                          %Set the heigh of the uicontrols.
    fontsize = 0.6*28.34*ui_h;                                              %Set the fontsize for all uicontrols.
    sp1 = 0.02*h;                                                           %Set the vertical spacing between axes and uicontrols.
    sp2 = 0.01*w;                                                           %Set the horizontal spacing between axes and uicontrols.
%     pos = [7*sp2,3*sp1,w-8*sp2,h-ui_h-5*sp1];                               %Set the position of the axes.
    pos = [7*sp2,7*sp1,w-8*sp2,h-ui_h-10*sp1];
    ax = axes('units','centimeters','position',pos,'box','on',...
        'linewidth',2);                                                     %Create axes for showing the log events histogram.
    obj = zeros(1,6);                                                       %Create a matrix to hold timescale uicontrol handles.
    str = {'Overall Hit Rate',...
        'Total Trial Count',...
        'Mean Peak Force',...
        'Median Peak Force',...
        'Trial Count',...
        'Hits in First 5 Minutes',...
        'Trials in First 5 Minutes',...
        'Max. Hits in Any 5 Minutes',...
        'Max. Trials in Any 5 Minutes',...
        'Max. Hit Rate in Any 5 Minutes',...
        'Min. Inter-Trial Interval (Smoothed)',...
        'Mean Peak Impulse',...
        'Median Peak Impulse',...
        'Peak Velocity',...
        'Latency to Hit'};                                             %List the available plots for the pull data.
    if any(strcmpi(devices{d},{'knob','lever'}))                            %If we're plotting knob data...
        str(2:3) = {'Mean Peak Angle','Median Peak Angle'};                 %Set the plots to show "angle" instead of "force".
    elseif ~any(strcmpi(devices{d},{'knob','lever','pull'}))                %Otherwise, if this isn't pull, knob, or lever data...
        str(2:3) = {'Mean Signal Peak','Median Signal Peak'};               %Set the plots to show "signal" instead of "peak force".
    end
    pos = [sp2, h-sp1-ui_h, 2*(w-6*sp2)/6, ui_h];                           %Set the position for the pop-up menu.
    obj(1) = uicontrol(fig,'style','popup','string',str,...
        'units','centimeters','position',pos,'fontsize',fontsize);      %Create pushbuttons for selecting the timescale.
    str = {'Session','Day','Week','Export'};                            %List the timescale labels.
    for i = 2:5                                                             %Step through the 3 timescales.
        pos = [i*sp2+i*(w-6*sp2)/6, h-sp1-ui_h, (w-6*sp2)/6, ui_h];         %Set the position for each pushbutton.
        obj(i) = uicontrol(fig,'style','pushbutton','string',str{i-1},...
            'units','centimeters','position',pos,'fontsize',fontsize);      %Create pushbuttons for selecting the timescale.
    end
    pos = [sp2, sp1, 2*(w-6*sp2)/6, ui_h];
    obj(6) = uicontrol(fig,'style','radiobutton','string','Concatenate Dates',...
        'units','centimeters','position',pos,'fontsize',fontsize);
    pos = [30*sp2, sp1, 2*(w-6*sp2)/6, ui_h];
    obj(7) = uicontrol(fig,'style','radiobutton','string','Group Animals',...
        'units','centimeters','position',pos,'fontsize',fontsize);
    set(obj(1),'callback',{@Set_Plot_Type,obj,data});                            %Set the callback for the pop-up menu.
    set(obj(2:4),'callback',{@Plot_Timeline,obj,[],data});                       %Set the callback for the timescale buttons.
    set(obj(5),'callback',{@Export_Data,ax,obj,data});                           %Set the callback for the export button.
    set(obj(6),'callback',{@TrainingDays,obj,data});
    set(obj(7),'callback',{@GroupAnimals,obj,data});
    set(fig,'userdata',data);                                            %Save the plot data to the figure's 'UserData' property.
    Plot_Timeline(obj(2),[],obj,[],data);                                        %Call the function to plot the session data in the figure.
    set(fig,'ResizeFcn',{@Resize,ax,obj});
end

%% This function is called when the user selects a plot type in the pop-up menu.
function Set_Plot_Type(~,~,obj,data)
i = strcmpi(get(obj,'fontweight'),'bold');                                  %Find the pushbutton with the bold fontweight.
Plot_Timeline(obj(i),[],obj,[],data);                                            %Call the subfunction to plot the data by the appropriate timeline.

%% This subfunction sorts the data into single-session values and sends it to the plot function.
function Plot_Timeline(hObject,~,obj,fid,data)
value = get(obj(6),'value');
GroupValue = get(obj(7),'value');
if GroupValue == 1 && isempty(get(obj(7),'userdata')) == 1;
    prompt = {'# of Groups', 'Group Names'};
    dlgtitle = 'Group Info';
    GroupInfo = inputdlg(prompt,dlgtitle, [1 20;1 50]);
    NumberGroups = str2double(GroupInfo{1});
    GroupTitle = strsplit(GroupInfo{2});
    rat_list = {data.rat};
    for i = 1:NumberGroups;
        Group_Selection = listdlg('PromptString',['Group: ' GroupTitle{i}],...
            'name','Subject Assignment',...
            'SelectionMode','multiple',...
            'listsize',[250 250],...
            'initialvalue',1:length(data),...
            'uh',25,...
            'ListString',rat_list);
        GroupData(i).name = GroupTitle(i);
        GroupData(i).rats = rat_list(Group_Selection);
    end
    set(obj(7),'userdata',GroupData);
end
if value == 1;
    Training_Days_Value = 'Training';
else
    Training_Days_Value = 'Dates';
end
TrialViewerData = data;
set(hObject,'fontweight','bold','foregroundcolor',[0 0.5 0]);               %Make this pushbutton's text bold.
set(setdiff(obj(2:4),hObject),'fontweight','normal','foregroundcolor','k'); %Make the other pushbutton's text normal and black.
fig = get(hObject,'parent');                                                %Grab the parent figure of the pushbutton.
data = get(fig,'userdata');                                                 %Grab the plot data from the figure's 'UserData' property.
i = find(hObject == obj);                                                   %Find the index of the button that called the function.
% if strcmpi(Training_Days_Value,'Dates') == 1;
    t = unique(horzcat(data.times));                                            %Horizontally concatenate all of the timestamps.
    if i == 2                                                                   %If the user wants to plot by session...
        t = [t; t + 0.00001]';                                                  %Set the time bounds to enclose only a single session.
    elseif i == 3                                                               %If the user wants to plot by day...
        t = unique(fix(t));                                                     %Find the unique truncated serial date numbers.
        t = [t; t + 1]';                                                        %Set the time bounds to go from the start of the day to the end of the day.
    else                                                                        %Otherwise, if the user wants to plot by week...
        t = [min(fix(t)), max(fix(t))];                                         %Find the first and last timestamp.
        i = find(strcmpi({'sun','mon','tue','wed','thu','fri','sat'},...
            datestr(t(1),'ddd')));                                              %Find the index for the day of the week of the first timestamp.
        t(1) = t(1) - i + 1;                                                    %Round down the timestamp to the nearest Sunday.
        t = t(1):7:t(2);                                                        %Find the timestamps for weekly spacing.
        t = [t; t + 7]';                                                        %Set the time bounds to go from the start of the week to the end of the week.
    end
% else
%     t = unique(horzcat(data.training_days));
%     if i == 2
%         t = [t; t+ 0.00001]';
%     end
% end
str = get(obj(1),'string');                                                 %Grab the strings from the pop-up menu.
i = get(obj(1),'value');                                                    %Grab the value of the pop-up menu.
str = str{i};                                                               %Grab the selected plot type.
plotdata = struct([]);                                                      %Create a structure to hold plot data.
for r = 1:length(data)                                                      %Step through each rat in the data structure.
    plotdata(r).rat = data(r).rat;                                          %Copy the rat name to the plot data structure.
    y = nan(1,size(t,1));                                                   %Pre-allocate a matrix to hold the data y-coordinates.
    s = cell(1,size(t,1));                                                  %Pre-allocate a cell array to hold the last stage of each time frame.
    n = cell(1,size(t,1));                                                  %Pre-allocate a cell array to hold the hit rate and trial count text.
    for i = 1:size(t,1)                                                     %Step through the specified time frames.
        j = data(r).times >= t(i,1) & data(r).times < t(i,2);               %Find all sessions within the time frame.
        if any(j)                                                           %If any sessions are found.
            if strcmpi(str,'overall hit rate')                              %If we're plotting overall hit rate...
                y(i) = nanmean(data(r).hitrate(j));                         %Grab the mean hit rate over this time frame.
            elseif strcmpi(str,'total trial count')                         %If we're plotting trial count...
                y(i) = nanmean(data(r).numtrials(j));                       %Grab the mean number of trials over this time frame.
            elseif any(strcmpi(str,{'median peak force',...
                    'median peak angle','median signal peak'}))             %If we're plotting the median signal peak...
                y(i) = nanmedian(data(r).peak(j));                          %Grab the mean signal peak over this time frame.
            elseif any(strcmpi(str,{'mean peak force',...
                    'mean peak angle','mean signal peak'}))                 %If we're plotting the mean signal peak...
                y(i) = nanmean(data(r).peak(j));                            %Grab the mean signal peak over this time frame.
            elseif strcmpi(str,'trial count')                               %If we're plotting number of trials....
                y(i) = nanmean(data(r).numtrials(j));                       %Grab the mean number of trials over this time frame.
            elseif strcmpi(str,'hits in first 5 minutes')                   %If we're plotting the hit count within the first 5 minutes.
                y(i) = nanmean(data(r).first_hit_five(j));                  %Grab the mean number of hits within the first 5 minutes over this time frame.
            elseif strcmpi(str,'trials in first 5 minutes')                 %If we're plotting the trial count within the first 5 minutes.
                y(i) = nanmean(data(r).first_trial_five(j));                %Grab the mean number of hits within the first 5 minutes over this time frame.
            elseif strcmpi(str,'max. hits in any 5 minutes')                %If we're plotting the maximum hit count within any 5 minutes.
                y(i) = nanmean(data(r).any_hit_five(j));                    %Grab the mean maximum number of hits within any 5 minutes over this time frame.
            elseif strcmpi(str,'max. trials in any 5 minutes')              %If we're plotting the maximum trial count within any 5 minutes.
                y(i) = nanmean(data(r).any_trial_five(j));                  %Grab the mean maximum number of trials within any 5 minutes over this time frame.
            elseif strcmpi(str,'max. hit rate in any 5 minutes')            %If we're plotting the maximum hit rate within any 5 minutes.
                y(i) = nanmean(data(r).any_hitrate_five(j));                %Grab the mean maximum hit rate within any 5 minutes over this time frame.
            elseif strcmpi(str,'min. inter-trial interval (smoothed)')      %If we're plotting the minimum inter-trial interval.
                y(i) = nanmean(data(r).min_iti(j));                         %Grab the mean minimum inter-trial interval over this time frame.
            elseif strcmpi(str,'median peak impulse')                       %If we're plotting the median signal impulse...
                y(i) = nanmedian(data(r).impulse(j));                       %Grab the mean signal impulse over this time frame.
            elseif strcmpi(str,'mean peak impulse')                         %If we're plotting the mean signal impulse...
                y(i) = nanmean(data(r).impulse(j));                         %Grab the mean signal impulse over this time frame.
            elseif strcmpi(str,'peak velocity');
                y(i) = nanmean(data(r).peak_velocity(j));
            elseif strcmpi(str,'latency to hit');
                y(i) = nanmean(data(r).latency(j));
            end
            temp = [nanmean(data(r).hitrate(j)),...
                nansum(data(r).numtrials(j))];                              %Grab the mean hit rate and total number of trials over this time frame.
            temp(1) = temp(1)*temp(2);                                      %Calculate the number of hits in the total number of trials.
            n{i} = sprintf('%1.0f hits/%1.0f trials',temp);                 %Create a string showing the number of hits and trials.
            j = find(j,1,'last');                                           %Find the last matching session.
            s{i} = data(r).stage{j};                                        %Save the last stage the rat ran on for this trime frame.            
        end
    end
    if strcmpi(Training_Days_Value,'Dates');
        plotdata(r).x = t(~isnan(y),:);                                         %Grab the daycodes at the start of each time frame.
    else
        plotdata(r).x = t(~isnan(y),:);
        plotdata(r).x_alt = 1:1:length(t(~isnan(y),:));
        plotdata(r).x_alt = [plotdata(r).x_alt; plotdata(r).x_alt + .2]';
    end      
    plotdata(r).y = y(~isnan(y))';                                          %Save only the non-NaN y-coordinates.
    plotdata(r).s = s(~isnan(y))';                                          %Save the stage information for each time frame.
    plotdata(r).n = n(~isnan(y))';                                          %Save the hit rate and trial information for each time frame.
end
ax = get(fig,'children');                                                   %Grab all children of the figure.
ax(~strcmpi(get(ax,'type'),'axes')) = [];                                   %Kick out all non-axes objects.
if isempty(fid)                                                             %If no text file handle was passed to this function...
    Make_Plot(plotdata,ax,str,TrialViewerData,obj);                                             %Call the subfunction to make the plot.
else                                                                        %Otherwise...
    [file, path] = uiputfile('*.xls','Save Spreadsheet');
    filename = [path file];
    for r = 1:length(plotdata)
        t = vertcat(plotdata(r).x);
        t = unique(t,'rows');
        t = cellstr(datestr(t(:,1)));
        meat = plotdata(r).y;
        name = cellstr(plotdata(r).rat);
        xlswrite(filename, {'Rat Name:'}, r,'A1');
        xlswrite(filename, name, r, 'B1');
        xlswrite(filename, {'Date/Time'},r,'A2');
        xlswrite(filename, {'Variable'},r,'B2');
        xlswrite(filename, t,r,'A3');
        xlswrite(filename, meat,r,'B3');
        clear t meat name
    end
%     t = vertcat(plotdata.x);                                                %Vertically concatenate all time-points.
%     t = unique(t,'rows');                                                   %Find all unique rows of the timepoints.
%     t = cellstr(datestr(t(:,1)));
%     meat = plotdata.y;
%     name = cellstr(plotdata.rat);
%     xlswrite(filename, {'Rat Name:'}, 1,'A1');
%     xlswrite(filename, name, 1, 'B1');
%     xlswrite(filename, {'Date/Time'},1,'A2');
%     xlswrite(filename, {'Variable'},1,'B2');
%     xlswrite(filename, t,1,'A3');
%     xlswrite(filename, meat,1,'B3');
    winopen(filename)
%     fprintf(fid,'%s,\t','DATE/TIME');                                       %Print a date column header.
%     for r = 1:length(plotdata)                                              %Step through the rats.
%         if r == length(plotdata)                                            %If this is the last rat...
%             fprintf(fid,'%s,\n',plotdata(r).rat);                           %Print the rat name followed by a carraige return.
%         else                                                                %Otherwise...
%             fprintf(fid,'%s,\t',plotdata(r).rat);                           %Print the rat name followed by a tab.
%         end
%     end
%     for i = 1:size(t,1)                                                     %Step through each time-point.
%         if rem(t(i,1),1) ~= 0                                               %If the timestamp is a fractional number of days...
%             fprintf(fid,'%s,\t',datestr(t(i,1),'mm/dd/yyyy - HH:MM'));      %Show the date and the time.
%         elseif t(i,2) - t(i,1) == 1                                         %If the timestamps only cover one day...
%             fprintf(fid,'%s,\t',datestr(t(i,1),'mm/dd/yyyy'));              %Show only the date.
%         else                                                                %Otherwise...
%             fprintf(fid,'%s,\t',[datestr(t(i,1),'mm/dd/yyyy') '-' ...
%                 datestr(t(i,2)-1,'mm/dd/yyyy')]);                           %Show the date range.
%         end
%         for r = 1:length(plotdata)                                          %Step through the rats.
%             j = plotdata(r).x(:,1) == t(i,1) & ...
%                 plotdata(r).x(:,1) == t(i,1);                               %Find any matching timepoint for this rat.
%             if any(j)                                                       %If any matching timepoint was found...
%                 fprintf(fid,'%1.3f,',plotdata(r).y(j));                     %Print the value for this date range.
%             else                                                            %Otherwise...
%                 fprintf(fid,'-,');                                          %Print a hyphen.
%             end
%             if r == length(plotdata)                                        %If this is the last rat...
%                 fprintf(fid,'\n');                                          %Print a carraige return.
%             else                                                            %Otherwise...
%                 fprintf(fid,'\t');                                          %Print a tab.
%             end
%         end
%     end
end

%% This subfunction sorts the data into daily values and sends it to the plot function.
function Export_Data(hObject,~,ax,obj,data)
output = questdlg(['Would you like to export the data as a spreadsheet '...
    'or figure image?'],'Data Type?','Spreadsheet','Image','Both','Both');  %Ask the user if they'd like to export the data as a spreadsheet or image.
fig = get(hObject,'parent');                                                %Grab the parent figure of the export button.
if any(strcmpi({'image','both'},output))                                    %If the user wants to save an image...
    [file, path] = uiputfile({'*.png';'*.jpg';'*.pdf'},...
        'Save Figure Image');                                               %Ask the user for a filename.
    if file(1) == 0                                                         %If the user clicked cancel...
        return                                                              %Skip execution of the rest of the function.
    end
    set(fig,'units','centimeters','resizefcn',[]);                          %Set the figure units to centimeters.
    pos = get(fig,'position');                                              %Get the figure position.
    temp = get(fig,'color');                                                %Grab the curret figure color.
    set(fig,'paperunits','centimeters','papersize',pos(3:4),...
        'inverthardcopy','off','paperposition',[0 0 pos(3:4)],'color','w'); %Set the paper size and paper position, in centimeters.
    w = pos(3);                                                             %Grab the width of the figure.
    h = pos(4);                                                             %Grab the height of the figure.
    ui_h = 0.07*h;                                                          %Set the height of all uicontrols.
    sp1 = 0.02*h;                                                           %Set the vertical spacing between axes and uicontrols.
    sp2 = 0.01*w;                                                           %Set the horizontal spacing between axes and uicontrols.
    pos = [7*sp2,3*sp1,w-8*sp2,h-4*sp1];                                    %Create an axes position matrix.
    
    set(ax,'units','centimeters','position',pos);                           %Expand the axes to fill the figure.
    set(obj,'visible','off');                                               %Make all of the other figures invisible.
    drawnow;                                                                %Immediately update the figure.    
    a = find(file == '.',1,'last');                                         %Find the start of the file extension.
    ext = file(a:end);                                                      %Grab the file extension.
    if strcmpi(ext,'.png')                                                  %If the user chose to save as a PNG...
        print(fig,[path file(1:a-1)],'-dpng');                              %Save the figure as a PNG file.
    elseif strcmpi(ext,'.jpg')                                              %If the user chose to save as a JPEG...
        print(fig,[path file(1:a-1)],'-djpeg');                             %Save the figure as a JPEG file.
    elseif strcmpi(ext,'.pdf')                                              %If the user chose to save as a PDF...
        print(fig,[path file(1:a-1)],'-dpdf');                              %Save the figure as a PDF file.
    end
%     pos = [7*sp2,3*sp1,w-8*sp2,h-ui_h-5*sp1];                               %Create an axes position matrix.
pos = [7*sp2,7*sp1,w-8*sp2,h-ui_h-10*sp1]; 
    set(ax,'units','centimeters','position',pos);                           %Reset the position of the axes.
    set(obj,'visible','on');                                                %Make all of the other figures visible again.
    i = strcmpi(get(obj,'fontweight'),'bold');                              %Find the pushbutton with the bold fontweight.
    Plot_Timeline(obj(i),[],obj,[],data);                                        %Call the subfunction to plot the data by the appropriate timeline.
    set(fig,'color',temp,'ResizeFcn',{@Resize,ax,obj});                     %Set the Resize function for the figure.
    drawnow;                                                                %Immediately update the figure.    
end
if any(strcmpi({'spreadsheet','both'},output))                              %If the user wants to save a spreadsheet...
%     temp = get(ax,'ylabel');                                                %Grab the handle for the axes y-label.
%     file = lower(get(temp,'string'));                                       %Grab the axes y-axis label.
%     file(file == ' ') = '_';                                                %Replace any spaces with underscores.
%     for i = '<>:"/\|?*().'                                                  %Step through all reserved characters.
%         file(file == i) = [];                                               %Kick out any reserved characters.
%     end
%     file = [file '_' datestr(now,'yyyymmdd')];                              %Add today's date to the default filename.
%     temp = lower(get(fig,'name'));                                          %Grab the figure name.
%     file = [temp(20:end) '_' file];                                         %Add the device name to the default filename.
%     [file, path] = uiputfile('*.csv','Save Spreadsheet',file);              %Ask the user for a filename.
%     if file(1) == 0                                                         %If the user clicked cancel...
%         return                                                              %Skip execution of the rest of the function.
%     end
%     fid = fopen([path file],'wt');                                          %Open a new text file to write the data to.   
    i = strcmpi(get(obj,'fontweight'),'bold');                              %Find the pushbutton with the bold fontweight.
fid = 1;
    Plot_Timeline(obj(i),[],obj,fid,data);                                       %Call the subfunction to write the data by the appropriate timeline.
%     fclose(fid);                                                            %Close the figure.
%     winopen([path file]);                                                   %Open the CSV file.
end

%% This section plots session/daily/weekly data in the specified axes.
function Make_Plot(plotdata,ax,str,TrialViewerData,obj)
Training_Days_Value = get(obj(6),'value');
GroupValue = get(obj(7),'value');
GroupData = get(obj(7),'userdata');
if Training_Days_Value == 0;
    Training_Days_Value = 'Dates';
end
fig = get(ax,'parent');                                                     %Grab the figure handle for the axes' parent.
set(fig,'units','centimeters');                                             %Set the figure's units to centimeters.
temp = get(fig,'position');                                                 %Grab the figure position.
h = temp(4);                                                                %Grab the figure height, in centimeters.
linewidth = 0.07*h;                                                          %Set the linewidth for the plots.
markersize = 0.75*h;                                                        %Set the marker size for the plots.
fontsize = 0.6*h;                                                           %Set the fontsize for all text objects.
cla(ax);                                                                    %Clear the axes.
if GroupValue == 0;
    colors = hsv(length(plotdata));                                             %Grab unique colors for all the rats.
else
    colors = hsv(length(GroupData));
end
hoverdata = struct([]);                                                     %Create an empty structure to hold data for the MouseHover function.
for r = 1:length(plotdata)                                                  %Step through each rat.
    if GroupValue == 1;
        q = 1; Counter = 1;
        while q == 1;
            compare = strcmpi(plotdata(r).rat,GroupData(Counter).rats);
            if any(compare) == 1;
                q = 0;
            else
                Counter = Counter + 1;
            end
        end
        c = Counter;
    else
        c = r;
    end
%     c = r;
    if strcmpi(Training_Days_Value, 'Dates') == 1;
        line(mean(plotdata(r).x,2),plotdata(r).y,'color',colors(c,:),...
            'linewidth',linewidth,'userdata',1,'parent',ax);                    %Show the rat's performance as a thick line.
    else
        line(mean(plotdata(r).x_alt,2),plotdata(r).y,'color',colors(c,:),...
            'linewidth',linewidth,'userdata',1,'parent',ax);                    %Show the rat's performance as a thick line.
    end
    for i = 1:size(plotdata(r).x,1)                                         %Step through each timepoint.
        if strcmpi(Training_Days_Value,'Dates') == 1;
            l(i) = line(mean(plotdata(r).x(i,:)),plotdata(r).y(i),...
                'markeredgecolor',colors(c,:),'linestyle','none',...
                'markerfacecolor',colors(c,:),'marker','.',...
                'linewidth',linewidth,'markersize',markersize,'userdata',2,...
                'parent',ax, 'ButtonDownFcn', {@TrialViewer,plotdata(r),i,TrialViewerData});                                                   %Mark each session with a unique marker.
        else
            l(i) = line(mean(plotdata(r).x_alt(i,:)),plotdata(r).y(i),...
                'markeredgecolor',colors(c,:),'linestyle','none',...
                'markerfacecolor',colors(c,:),'marker','.',...
                'linewidth',linewidth,'markersize',markersize,'userdata',2,...
                'parent',ax, 'ButtonDownFcn', {@TrialViewer,plotdata(r),i,TrialViewerData});
        end
        if strcmpi(Training_Days_Value,'Dates') == 1;
            hoverdata(end+1).xy = [mean(plotdata(r).x(i,:)),plotdata(r).y(i)];  %Save the x- and y-coordinates.
        else
            hoverdata(end+1).xy = [mean(plotdata(r).x_alt(i,:)),plotdata(r).y(i)];  %Save the x- and y-coordinates.
        end
        if rem(plotdata(r).x(i,1),1) ~= 0                                   %If the timestamp is a fractional number of days...
            temp = datestr(plotdata(r).x(i,1),'mm/dd/yyyy, HH:MM');         %Show the date and the time.
        elseif plotdata(r).x(i,2) - plotdata(r).x(i,1) == 1                 %If the timestamps only cover one day...
            temp = datestr(plotdata(r).x(i,1),'mm/dd/yyyy');                %Show only the date.
        else                                                                %Otherwise...
            temp = [datestr(plotdata(r).x(i,1),'mm/dd/yyyy') '-' ...
                datestr(plotdata(r).x(i,2)-1,'mm/dd/yyyy')];                %Show the date range.
        end
        hoverdata(end).txt = {plotdata(r).rat,plotdata(r).s{i},...
            temp,plotdata(r).n{i}};                                         %Save the rat's name, stage, date, and hit rate/num trials.
        hoverdata(end).handles = l;                                         %Save the line handle for the point.
    end
end
if strcmpi(Training_Days_Value,'Dates') == 1;
    temp = vertcat(hoverdata.xy);                                               %Grab all of the datapoint coordinates.
    x = [min(temp(:,1)), max(temp(:,1))];                                       %Find the minimim and maximum x-values.
    x = x + [-0.05,0.05]*(x(2) - x(1));                                         %Add some padding to the x-axis limits.
else
    temp = vertcat(plotdata.x_alt);
    x = [min(temp(:,1)), max(temp(:,1))];
    x = x + [-.05,.05]*(x(2) - x(1));
end
if length(x) < 2 || any(isnan(x))                                           %If there are any missing x values.
    x = now + [-1,1];                                                       %Set arbitrary x-axis limits.
end
if x(1) == x(2)                                                             %If the x-axis limits are the same...
    x = x + [-1,1];                                                         %Add a day to each side of the single point.
end
xlim(ax,x);                                                                 %Set the x-axis limits.
if strcmpi(Training_Days_Value,'Dates') == 1;
    y = [min(temp(:,2)), max(temp(:,2))];                                       %Find the minimim and maximum x-values.
    y = y + [-0.05,0.05]*(y(2) - y(1));                                         %Add some padding to the y-axis limits.
else
    temp = vertcat(plotdata.y);
    y = [min(temp(:,1)), max(temp(:,1))];
    y = y + [-0.05,0.05]*(y(2) - y(1));
end
if length(y) < 2 || any(isnan(y))                                           %If there are any missing y values.
    y = [0,1];                                                              %Set arbitrary y-axis limits.
end
if y(1) == y(2)                                                             %If the y-axis limits are the same...
    y = y + [-0.1, 0.1];                                                    %Add 0.1 to each side of the single point.
end
ylim(ax,y);                                                                 %Set the y-axis limits
if strcmpi(Training_Days_Value,'Dates') == 1;
    XTickLabelTemp = datestr(get(ax,'xtick'),'mm/dd');
else
    XTickLabelTemp = get(ax,'xtick');
end
set(ax,'xticklabel',XTickLabelTemp,...
    'fontsize',fontsize,'fontweight','bold','linewidth',linewidth);         %Show the date for each x-axis tick.
ylabel(ax,str,'fontweight','bold','fontsize',1.1*fontsize);                 %Label the x-axis.
if GroupValue == 1;
    for i = 1:length(GroupData);
        legendnames(i) = GroupData(i).name;
    end
    [legh,objh,outh,outm] = legend(legendnames,'Location','southeast');
%     lineh = findobj(objh,'type','line');
    texth = findobj(objh,'type','text');
    set(outh,'Linestyle','-','Marker','none','Linewidth',1);
    set(texth,'Fontsize',10, 'Fontweight', 'normal');
    set(legh,'Fontsize', 14);
%     lineh = findobj(lineh,'Marker','none');
    
else
    legend('hide')
end
temp = get(ax,'ytick');                                                     %Grab the y-axis ticks.
for i = 1:length(temp)                                                      %Step through the y-axis ticks.
    temp(i) = line(xlim(ax),temp(i)*[1,1],'color',[0.75 0.75 0.75],...
        'linestyle','--','linewidth',0.5*linewidth,'userdata',3,...
        'parent',ax);                                                       %Draw a gridline at each y-tick.
end
uistack(temp,'bottom');                                                     %Send all grid lines to the bottom.
txt = text(x(1),y(1),' ','fontsize',fontsize,'margin',2,...
    'backgroundcolor','w','edgecolor','k','visible','off',...
    'verticalalignment','bottom','horizontalalignment','center',...
    'userdata',4);                                                          %Create a text object for labeling points.
set(fig,'WindowButtonMotionFcn',{@MouseHover,ax,hoverdata,txt});            %Set the mouse hover function for the figure.
% set(fig,'ButtonDownFcn', {@TrialViewer,ax,hoverdata,txt});

%% This function executes while the mouse hovers over the axes of an interactive figure.
function MouseHover(~,~,ax,data,txt)
xy = get(ax,'CurrentPoint');                                                %Grab the current mouse position in the axes.
xy = xy(1,1:2);                                                             %Pare down the x-y coordinates matrix.
a = [get(ax,'xlim'), get(ax,'ylim')];                                       %Grab the x- and y-axis limits.
if xy(1) >= a(1) && xy(1) <= a(2) && xy(2) >= a(3) && xy(2) <= a(4)         %If the mouse was clicked inside the axes...
    fig = get(ax,'parent');                                                 %Grab the parent figure of the axes.
    set(fig,'units','centimeters');                                         %Set the figure units to centimeters
    pos = get(fig,'position');                                              %Grab the current figure size, in centimeters.
    markersize = 0.75*pos(4);                                               %Set the marker size for the plots.
    xy = (xy - a([1,3]))./[a(2) - a(1), a(4) - a(3)];                       %Normalize the mouse x-y coordinates.
    temp = vertcat(data.xy);                                                %Vertically concatenate the point x-y coordinates.
    temp(:,1) = (temp(:,1) - a(1))/(a(2) - a(1));                           %Normalize the point x coordinates.
    temp(:,2) = (temp(:,2) - a(3))/(a(4) - a(3));                           %Normalize the point 3 coordinates.
    set(ax,'units','centimeters');                                          %Set the axes position units to centimeters.
    a = get(ax,'position');                                                 %Get the axes position.
    xy = xy.*a(3:4);                                                        %Calculate the axes position in centimeters.
    temp(:,1) = a(3)*temp(:,1);                                             %Calculate the point x coordinates in centimeters
    temp(:,2) = a(4)*temp(:,2);                                             %Calculate the point y coordinates in centimeters
    d = [xy(1) - temp(:,1), xy(2) - temp(:,2)];                             %Calculate the x-y distances from the mouse to the points.
    d = sqrt(sum(d.^2,2));                                                  %Find the straight-line distance to each point.
    if any(d <= 0.5)                                                        %If we're within half a centimeter of any data point...
        i = find(d == min(d),1,'first');                                    %Find the closest point.
        xy = data(i).xy;                                                    %Grab the data point x-y coordinate.
        temp = xlim(ax);                                                    %Grab the x-axis limits.
        if xy(1) < temp(1) + 0.25*(temp(2) - temp(1))                       %If the data point is in the left-most quartile...
            set(txt,'horizontalalignment','left');                          %Set the horizontal alignment to left-hand.
        elseif xy(1) > temp(1) + 0.75*(temp(2) - temp(1))                   %If the data point is in the right-most quartile...
            set(txt,'horizontalalignment','right');                         %Set the horizontal alignment to right-hand.
        else                                                                %Otherwise...
            set(txt,'horizontalalignment','center');                        %Set the horizontal alignment to centered.
        end
        if xy(2) > mean(ylim(ax))                                           %If the data point is in the top half...
            xy(2) = xy(2) - 0.05*range(ylim);                               %Adjust the y coordinate to place the text below the data point.
            set(txt,'verticalalignment','top');                             %Set the vertical alignment to top.
        else                                                                %Otherwise...
            xy(2) = xy(2) + 0.05*range(ylim);                               %Adjust the y coordinate to place the text above the data point.
            set(txt,'verticalalignment','bottom');                          %Set the vertical alignment to bottom.
        end
        temp = get(txt,'position');                                         %Grab the current text object position.
        str = get(txt,'string');                                            %Grab the rat's name.
        if ~isequal(xy,temp) || ~strcmpi(str,data(i).txt)                   %If the current label is incorrect...
            set(txt,'position',xy,'string',data(i).txt,'visible','on');     %Update the position, string, and visibility of the text object.
            set(data(i).handles,'markersize',2*markersize);                 %Make the selected marker larger.
            set(setdiff([data.handles],data(i).handles),...
                'markersize',markersize);                                   %Make the other markers smaller.
        end
    else                                                                    %Otherwise...
        set([data.handles],'markersize',markersize);                        %Make all markers smaller.
        set(txt,'visible','off');                                           %Make the text object invisible.
    end
end

%% This function is called whenever the main figure is resized.
function Resize(hObject,~,ax,obj)
set(hObject,'units','centimeters');                                         %Set the figure units to centimeters
pos = get(hObject,'position');                                              %Grab the current figure size, in centimeters.
w = pos(3);                                                                 %Grab the width of the figure.
h = pos(4);                                                                 %Grab the height of the figure.
ui_h = 0.07*h;                                                              %Set the heigh of all uicontrols.
sp1 = 0.02*h;                                                               %Set the vertical spacing between axes and uicontrols.
sp2 = 0.01*w;                                                               %Set the horizontal spacing between axes and uicontrols.
fontsize = 0.6*28.34*ui_h;                                                  %Set the fontsize for all uicontrols.
% pos = [7*sp2,3*sp1,w-8*sp2,h-ui_h-5*sp1];   
pos = [7*sp2,7*sp1,w-8*sp2,h-ui_h-10*sp1];
%Create an axes position matrix.
set(ax,'units','centimeters','position', pos);                              %Reset the position of the axes.
pos = [sp2, h-sp1-ui_h, 2*(w-6*sp2)/6, ui_h];                               %Create the position matrix for the pop-up menu.
set(obj(1),'units','centimeters','position',pos,'fontsize',fontsize);       %Set the position of the pop-up menu.
for i = 2:5                                                                 %Step through the timescale and export pushbuttons.      
    pos = [i*sp2+i*(w-6*sp2)/6, h-sp1-ui_h, (w-6*sp2)/6, ui_h];             %Create the position matrix for each pushbutton.
    set(obj(i),'units','centimeters','position',pos,'fontsize',fontsize);   %Set the pushbutton position and fontsize.
end
pos = [sp2, sp1, 2*(w-6*sp2)/6, ui_h];
set(obj(6),'units','centimeters','position',pos,'fontsize',fontsize);
pos = [30*sp2, sp1, 2*(w-6*sp2)/6, ui_h];
set(obj(7),'units','centimeters','position',pos,'fontsize',fontsize);
linewidth = 0.07*h;                                                          %Set the linewidth for the plots.
markersize = 0.75*h;                                                        %Set the marker size for the plots.
fontsize = 0.6*h;                                                           %Set the fontsize for all text objects.
obj = get(ax,'children');                                                   %Grab all children of the axes.
temp = vertcat(get(obj,'userdata'));                                        %Vertically concatenate the 'UserData' from each objects.
temp = vertcat(temp{:});                                                    %Change the 'UserData' from a cell array to matrix.
i = strcmpi(get(obj,'type'),'line') & temp == 1;                            %Find all plot line objects.
set(obj(i),'linewidth',linewidth);                                          %Set the linewidth for all plot line objects.
i = strcmpi(get(obj,'type'),'line') & temp == 2;                            %Find all plot marker objects.
set(obj(i),'markersize',markersize);                                        %Set the markersize for all plot line objects.
i = strcmpi(get(obj,'type'),'line') & temp == 3;                            %Find all grid line objects.
set(obj(i),'linewidth',0.5*linewidth);                                      %Set the linewidth for all grid line objects.
i = strcmpi(get(obj,'type'),'text');                                        %Find all text objects.
set(obj(i),'fontsize',fontsize);                                            %Set the font size for all text objects.
set(ax,'fontsize',fontsize,'linewidth',linewidth);                          %Set the axes linewidth and fontsize.
temp = get(ax,'ylabel');                                                    %Grab the y-axis label handle for the axes.
set(temp,'fontsize',1.1*fontsize);                                          %Set the font size for y-axis label.

%% This function is called whenever a point is selected in General Analysis
function TrialViewer(hObject,~,plotdata,i,TrialViewerData)
choice = questdlg('Would you like to view each trial of this session?',...
    'Analyze Data',...
    'Yes', 'No', 'No');
switch choice
    case 'Yes'
        Files = TrialViewerData.files;
        a = find(Files{i} == '\',1,'last');                                     %Find the last forward slash in the filename.
        files = Files{i}(a+1:end);                                               %Grab the filename minus the path.
        path = Files{i}(1:a);
        MotoTrak_Supination_Viewer_Edit(files,path)
end

%% This subfunction finds peaks in the signal, accounting for equality of contiguous samples.
function [pks,i] = PeakFinder(signal,smoothsize,init)
temp = boxsmooth(signal,smoothsize);                                        %Smooth the signal.
i = find(temp(2:end) - temp(1:end-1) > 0) + 1;                              %Find each point that's greater than the preceding point.
j = find(temp(1:end-1) - temp(2:end) >= 0);                                 %Find each point that's greater than or equal to the following point.
i = intersect(i,j);                                                         %Find any points that meet both criteria.
j = (signal(i) < init);                                                     %Find all peaks smaller than the initiation threshold.
i(j) = [];                                                                  %Kick out all sub-initiation peaks.
pks = signal(i);                                                            %Grab the value of the signal at each peak.