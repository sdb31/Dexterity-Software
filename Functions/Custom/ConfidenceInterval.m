function ConfidenceInterval(varargin)

[temp_pre_files pre_path] = uigetfile({'*.MotoTrak;*.ArdyMotor'},...
    'Select Animal(s) Data', ...
    'multiselect','on');                                         %Have the user pick an input *.ArdyMotor file or pre_files.
if ~iscell(temp_pre_files) && temp_pre_files(1) == 0                              %If no file was selected...
    return                                                              %Exit the function.
end
cd(pre_path);                                                           %Change the current directory to the specified folder.
if ischar(temp_pre_files)                                                    %If only one file was selected...
    temp_pre_files = {temp_pre_files};                                            %Convert the string to a cell array.
end
try                                                                     %Try to read in the data file...
    [~,filename,ext] = fileparts(temp_pre_files{1});
    switch ext
        case '.ArdyMotor'
            temp = ArdyMotorFileRead(temp_pre_files{1});
        case '.MotoTrak'
            temp = MotoTrakFileRead(temp_pre_files{1});
            temp = MotoTrak_to_ArdyMotor(temp);
            temp.rat = temp.subject;
    end
    %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
catch err                                                               %If an error occurs...
    warning(['ERROR READING: ' files{1}]);                              %Show which file had a read problem...
    warning(err.message);                                               %Show the actual error message.
end
% temp = ArdyMotorFileRead(temp_pre_files{1});
animal_pre = temp.rat;

%This code sorts the files by daycode in ardymotor function
num_sessions = length(temp_pre_files);
for f = 1:num_sessions
    try                                                                     %Try to read in the data file...
        [~,filename,ext] = fileparts(temp_pre_files{f});
        switch ext
            case '.ArdyMotor'
                data = ArdyMotorFileRead(temp_pre_files{f});
                pre_daycodes(f) = data.daycode;
            case '.MotoTrak'
                data = MotoTrakFileRead(temp_pre_files{f});
                data = MotoTrak_to_ArdyMotor(data);
                data.rat = data.subject;
                pre_daycodes(f) = data.start_time;
        end
        %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
    catch err                                                               %If an error occurs...
        warning(['ERROR READING: ' files{1}]);                              %Show which file had a read problem...
        warning(err.message);                                               %Show the actual error message.
    end
%     data = ArdyMotorFileRead(temp_pre_files{f});
%     pre_daycodes(f) = data.daycode;
end
[B,IX] = sort(pre_daycodes);
for f = 1:num_sessions
    pre_files(f) = temp_pre_files(IX(f));
end

try
    handles.(animal_pre).pre.files = pre_files;                             %Store the pre_files
    handles.(animal_pre).pre.path = pre_path;
catch err
    animal_pre = 'TempRat';
    
    handles.(animal_pre).pre.files = pre_files;                             %Store the pre_files
    handles.(animal_pre).pre.path = pre_path;
end


handles.animals = animal_pre;


animal = handles.animals;

%% Do pre-lesion file processing
% handles.(animal).pre.data = TransformKnobData( ...
%     handles.(animal).pre.files, handles.(animal).pre.path);
handles.(animal).pre.trialdata = TransformTrialData( ...
    handles.(animal).pre.files, handles.(animal).pre.path);

%% Set up axes for plotting
handles.analysis = 'data';

handles.(animal).pre.session_select = 1;
handles.(animal).post.session_select = 1;

PlotTrial(handles, 'pre');

function PlotTrial(handles, preorpost)

% cla(handles.(preorpost).trial_axes);
figure; clf;

trial_start = 100;
trial_end = 300;

s_pre = handles.(handles.animals).pre.session_select;
% s_post = handles.(animal).post.session_select;


handles.pre.max_signal = max(handles.(handles.animals).pre.trialdata.mean_plot(s_pre,:));
% handles.post.max_signal = max(handles.(animal).post.trialdata.mean_plot(s_post,:));
% max_signal = max(handles.pre.max_signal, handles.post.max_signal);
max_signal = max(handles.pre.max_signal);

y_max = 1.2*max_signal;


y_min = -(0.025*y_max);
s = handles.(handles.animals).(preorpost).session_select;

%Plot confidence interval fill
% fill([1:500,500:-1:1], ... 
%     [handles.(handles.cur_animal).(preorpost).trialdata.upper(s,:), ... 
%     fliplr(handles.(handles.cur_animal).(preorpost).trialdata.lower(s,:))], ... 
%     [0 0.5 0], 'parent', handles.(preorpost).trial_axes);
fill([1:500,500:-1:1], ... 
    [handles.(handles.animals).(preorpost).trialdata.upper(s,:), ... 
    fliplr(handles.(handles.animals).(preorpost).trialdata.lower(s,:))], ... 
    [0 0.5 0]);
hold on;

%plot mean data line
% line((1:500),handles.(handles.cur_animal).(preorpost).trialdata.mean_plot(s,:), ... 
%     'Color', 'k', 'linewidth', 2, 'parent', handles.(preorpost).trial_axes);
line((1:500),handles.(handles.animals).(preorpost).trialdata.mean_plot(s,:), ... 
    'Color', 'k', 'linewidth', 2);
hold on;

%Plot max pre signal line
line([0 500], [handles.pre.max_signal handles.pre.max_signal], ...
    'linestyle', '--', 'color', 'b');
%Label max pre signal line
text(450, handles.pre.max_signal, ['Max Signal: ' num2str(handles.pre.max_signal)], ...
    'horizontalalignment','center','verticalalignment','middle', ...
    'margin',2,'edgecolor','k',...
    'backgroundcolor','w', 'fontsize',8, 'FontWeight', 'bold');

% %Plot max post signal line
% line([0 500], [handles.post.max_signal handles.post.max_signal], ...
%     'linestyle', '--', 'color', [1 0.4 0.7]);
% %Label max post signal line
% text(450, handles.post.max_signal, ['Post Max: ' num2str(handles.post.max_signal)],...
%     'horizontalalignment','center','verticalalignment','middle', ...
%     'margin',2,'edgecolor','k',...
%     'backgroundcolor','w', 'fontsize',8,'FontWeight', 'bold');


%Plot x axis line 
line([0 500], [0 0], 'linewidth', 2, 'color', 'k');

%plot hit window start line
line([trial_start trial_start], [0 y_max], ...
    'Color', 'k', 'linestyle', '--', 'linewidth', 1);

%Plot hit window end line
line([trial_end trial_end], [0 y_max], ...
    'Color', 'k', 'linestyle', '--', 'linewidth', 1);

%plot Hit Window line up top
line([trial_start trial_end], [0.95*y_max 0.95*y_max], ...
    'Color', 'k', 'linestyle', '--', 'linewidth', 1);

%Plot hit window text box
text(trial_end - trial_start,0.95*y_max, ... 
    'Hit Window','margin',2,'edgecolor','k',...
    'backgroundcolor','w', 'fontsize',10,...
    'horizontalalignment','center','verticalalignment','middle', ...
    'FontWeight', 'bold'); 

%Label axes

% if (strcmpi(preorpost, 'pre'))
%     title_label = 'Pre';
% else
%     title_label = 'Post';
% end

title([handles.animals ': Mean Trial with 95% Confidence Interval'], ...
    'FontSize', 12, 'FontWeight', 'bold');
ylabel('Degrees', 'FontWeight', 'bold', ...
    'FontSize', 12);
xlabel('Time (10ms)', 'FontWeight', 'bold');


%% This function is called whenever the main figure is resized
function Resize(hObject,~)
handles = guidata(hObject);                                                 %Grab the handles structure from the GUI.
pos = get(handles.fig,'position');                                          %Grab the main figure position.
% ylabel('Test','parent',handles.pre.axes,'fontsize',0.75*pos(4),...
%     'rotation',0,'verticalalignment','middle',...
%     'horizontalalignment','right');                                    
% ylabel('TestPost','parent',handles.post.axes,'fontsize',0.75*pos(4),...
%     'rotation',0,'verticalalignment','middle',...
%     'horizontalalignment','right');     

set([handles.label,handles.savebutton],'fontsize',0.75*pos(4));             %Update the trial label and savebutton fontsize.

objs_pre = get(handles.pre.axes,'children');                                  %Grab all children of the force axes.
objs_post = get(handles.post.axes,'children');

% objs_pre(~strcmpi('text',get(objs_pre,'type'))) = [];                               %Kick out all non-text objects.
% objs_post(~strcmpi('text',get(objs_post,'type'))) = [];                               %Kick out all non-text objects.


set(objs_pre,'fontsize',0.5*pos(4));                                            %Update the fontsize of all text objects.
set(objs_post, 'fontsize', 0.5*post(4));

ylabel('TestPre)','parent',handles.pre.axes,'fontsize',0.75*pos(4));     %Label the force signal.
xlabel('Trial','parent',handles.pre.axes,'fontsize',0.75*pos(4));     %Label the time axis.

ylabel('TestPost','parent',handles.post.axes,'fontsize',0.75*pos(4));     %Label the force signal.
xlabel('Trial','parent',handles.post.axes,'fontsize',0.75*pos(4));     %Label the time axis.

%% This subfunction finds peaks in the signal, accounting for equality of contiguous samples.
function [pks, sig] = Knob_Peak_Finder(signal)
    %This code finds and kicks out peaks that have a std dev between 
    %them less than 1
    
    smoothed_signal = boxsmooth(signal);                                        %smooth out the trial signal
    [pks, sig] = findpeaks(smoothed_signal, 'MINPEAKHEIGHT', 5, ...
        'MINPEAKDISTANCE', 10);                                            %Find local maximma
    n = length(pks);
    j = 1;
    if n>1
        while j <= n-1
            if (abs(pks(j)-pks(j+1)) <= 5)                                 % if the diff between 2 peaks is less than or equal to 5
                start_sig = sig(j);
                end_sig = sig(j+1);

                signal_interest = smoothed_signal(start_sig:end_sig);
                deviation_signal = std(signal_interest);

                if deviation_signal < 1
                    pks(j+1) = [];
                    sig(j+1) = [];
                    j = j-1;
                end

            end
            n = length(pks);
            j = j+1;
        end
    end


%% This subfunction returns the daycode (1-365) for a given date.
function d = daycode(date)
date = datevec(date);                                                       %Convert the serial date number to a date vector.
year = date(1);                                                             %Pull the year out of the date vector.
month = date(2);                                                            %Pull out the month.
day = date(3);                                                              %Pull out the day.
if year/4 == fix(year/4);                                                   %If the year is a leap year...
    numDays = [31 29 31 30 31 30 31 31 30 31 30 31];                        %Include 29 days in February.
else                                                                        %Otherwise...
	numDays = [31 28 31 30 31 30 31 31 30 31 30 31];                        %Include 28 days in February.
end
date = sum(numDays(1:(month-1)));                                           %Sum the days in the preceding months...
d = date + day;                                                             %...and add the day of the specified month.

function X = boxsmooth(X,wsize)
%Box smoothing function for 2-D matrices.

%X = BOXSMOOTH(X,WSIZE) performs a box-type smoothing function on 2-D
%matrices with window width and height equal to WSIZE.  If WSIZE isn't
%given, the function uses a default value of 5.

if (nargin < 2)                                                             %If the use didn't specify a box size...
    wsize = 5;                                                              %Set the default box size to a 5x5 square.
end     
if (nargin < 1)                                                             %If the user entered no input arguments...
   error('BoxSmooth requires 2-D matrix input.');                           %Show an error.
end

if length(wsize) == 1                                                       %If the user only inputted one dimension...
    rb = round(wsize);                                                      %Round the number of row bins to the nearest integer.
    cb = rb;                                                                %Set the number of column bins equal to the number of row bins.
elseif length(wsize) == 2                                                   %If the user inputted two dimensions...
    rb = round(wsize(1));                                                   %Round the number of row bins to the nearest integer.
    cb = round(wsize(2));                                                   %Round the number of column bins to the nearest integer.
else                                                                        %Otherwise, if the 
    error('The input box size for the boxsmooth can only be a one- or two-element matrix.');
end

w = ones(rb,cb);                                                            %Make a matrix to hold bin weights.
if rem(rb,2) == 0                                                           %If the number of row bins is an even number.
    rb = rb + 1;                                                            %Add an extra bin to the number of row bins.
    w([1,end+1],:) = 0.5;                                                   %Set the tail bins to have half-weight.
end
if rem(cb,2) == 0                                                           %If the number of column bins is an even number.
    cb = cb + 1;                                                            %Add an extra bin to the number of row bins.
    w(:,end+1) = w(:,1);                                                    %Make a new column of weights with the weight of the first column.
    w(:,[1,end]) = 0.5*w(:,[1,end]);                                        %Set the tail bins to have half-weight.
end

[r,c] = size(X);                                                            %Find the number of rows and columns in the input matrix.
S = nan(r+rb-1,c+cb-1);                                                     %Pre-allocate an over-sized matrix to hold the original data.
S((1:r)+(rb-1)/2,(1:c)+(cb-1)/2) = X;                                       %Copy the original matrix to the center of the over-sized matrix.

temp = zeros(size(w));                                                      %Pre-allocate a temporary matrix to hold the box values.
for i = 1:r                                                                 %Step through each row of the original matrix.
    for j = 1:c                                                             %Step through each column of the original matrix.
        temp(:) = S(i:(i+rb-1),j:(j+cb-1));                                 %Pull all of the bin values into a temporary matrix.
        k = ~isnan(temp(:));                                                %Find all the non-NaN bins.
        X(i,j) = sum(w(k).*temp(k))/sum(w(k));                              %Find the weighted mean of the box and save it to the original matrix.
    end
end

%% This function analyzes the data 
function [knob_data] = TransformKnobData(files, path)

parameters = {'mean_distance', 'latency_to_hit', 'peak_velocity', ...
    'peak_acceleration', 'num_turns', 'mean_turn_distance'};


%% Find the maximum length across all files
max_session = 0;
cd(path);

num_sessions = length(files);
knob_data.num_sessions = num_sessions;

for f = 1:knob_data.num_sessions
    data = ArdyMotorFileRead(files{f});
    temp = length(data.trial);
    if temp > max_session
        max_session = temp;
    end
end



%% Declare all variables

for p = 1:length(parameters)
    % Create matricies to store the totals trial variables
    knob_data.(parameters{p}).sessions = nan(num_sessions, max_session);
    %Create matricies to store the means of individual sessions
    knob_data.(parameters{p}).sessionmeans = nan(num_sessions, 1);
    %Create matricies to store the combined session data
    knob_data.(parameters{p}).sessionscombined = nan(1,num_sessions);
end

knob_data.threshold.combined = nan(1, num_sessions);

for f = 1:num_sessions                                                     %Step through each file.
    
    cd(path);
    data = ArdyMotorFileRead(files{f});                                     %Read in the data from the *.ArdyMotor file.
    if ~isfield(data,'trial') || isempty(data.trial)                        %If there's no trials in this data file...
        warning('ARDYMOTOR2TEXT:NoTrials',['WARNING FROM '...
            'ARDYMOTOR2TEXT: The file "' file{f} '" has zero trials '...
            'and will be skipped.']);                                       %Show a warning.
        continue                                                            %Skip to the next file.
    end
    
    knob_data.sessions(f,1) = length(data.trial);

    for t = 1:length(data.trial)                                            %Step through each trial.
        a = find((data.trial(t).sample_times >= 0 & ...
        data.trial(t).sample_times < 1000*data.trial(t).hitwin));           %Find all samples within the hit window.
    
        signal = data.trial(t).signal(a);                                   %Define the signal as the raw signal
        
        data.trial(t).range = range(signal);                                %Find the hit window range of each trial signal.
        data.trial(t).max = max(signal);                                    %Find the hit window maximum of each trial signal.
        data.trial(t).min = min(signal);                                    %Find the hit window minimum of each trial signal.
        
        [pks, sig] = Knob_Peak_Finder(signal);                              %Find the peaks in the signal
        data.trial(t).num_turns = length(sig);                              %Calculate number of turns
        data.trial(t).mean_turn_distance = mean(pks);                       %Take the mean peak distance to find mean turn distance in trial
        
        if data.trial(t).num_turns == 0                                     %If a turn was initiated in the previous session, and not registered
            data.trial(t).num_turns = 1;                                    %in this session, then count it as 1
        end
         
        smooth_trial_signal = boxsmooth(signal);                            %Smooth the signal 
        
        smooth_knob_velocity = boxsmooth(diff(smooth_trial_signal));        %Boxsmooth the velocity signal
        data.trial(t).peak_velocity = max(smooth_knob_velocity);            %Calculate the max of the velocity signal
        
        knob_acceleration = boxsmooth(diff(smooth_knob_velocity));          %Calculate the acceleration by differentiating the velocity signal
        data.trial(t).peak_acceleration = max(knob_acceleration);           %Find the max of the acceleration signal
        
        if (data.trial(t).outcome == 72)                                    %If it was a hit
            hit_time = find(data.trial(t).signal >= ...                     %Calculate the hit time
                data.trial(t).thresh,1);                                   
            data.trial(t).latency_to_hit = hit_time;                        %Set this hit_time to the latency_to_hit structure variable
        else
            data.trial(t).latency_to_hit = NaN;                             %If trial resulted in a miss, then set latency to hit to NaN
        end
        
        %Put each parameter of individual trial in its appropriate array,
        %store each trial from a session in consecutive columns,
        %and each session in consecutive rows
        
        
            knob_data.mean_distance.sessions(f,t) = data.trial(t).max;                   
            knob_data.threshold.sessions(f,t) = data.trial(t).thresh;
            knob_data.latency_to_hit.sessions(f,t) = data.trial(t).latency_to_hit;
            knob_data.peak_velocity.sessions(f,t) = data.trial(t).peak_velocity;
            knob_data.peak_acceleration.sessions(f,t) = data.trial(t).peak_acceleration;
            knob_data.num_turns.sessions(f,t) = data.trial(t).num_turns;
            knob_data.mean_turn_distance.sessions(f,t) = data.trial(t).mean_turn_distance;
        
    end
    
    
    
    %Store the total mean of the session in consecutive rows
    knob_data.mean_distance.sessionmeans(f,1) = nanmean(knob_data.mean_distance.sessions(f,:));
    knob_data.latency_to_hit.sessionmeans(f,1) = nanmean(knob_data.latency_to_hit.sessions(f,:));
    knob_data.peak_velocity.sessionmeans(f,1) = nanmean(knob_data.peak_velocity.sessions(f,:));
    knob_data.peak_acceleration.sessionmeans(f,1) = nanmean(knob_data.peak_acceleration.sessions(f,:));
    knob_data.num_turns.sessionmeans(f,1) = nanmean(knob_data.num_turns.sessions(f,:));
    knob_data.mean_turn_distance.sessionmeans(f,1) = nanmean(knob_data.mean_turn_distance.sessions(f,:));

end

for p = 1:length(parameters)
    prev_sessions = 1;
    
    %Concatenate all trials from consecutive sessions in to one long array
    for s = 1:num_sessions
        knob_data.(parameters{p}).sessionscombined(prev_sessions:(knob_data.sessions(s)+prev_sessions-1)) ...       %Throw everything in to one array
            = knob_data.(parameters{p}).sessions(s,1:knob_data.sessions(s));
        if (strcmpi((parameters{p}), 'mean_distance') == 1)
            knob_data.threshold.sessionscombined(prev_sessions:(knob_data.sessions(s)+prev_sessions-1)) ...
                = knob_data.threshold.sessions(s,1:knob_data.sessions(s));
        end
        prev_sessions = knob_data.sessions(s) + prev_sessions;
    end
    
    %Calculate the mean of concatenated matrix
    knob_data.(parameters{p}).combinedmean = nanmean(knob_data.(parameters{p}).sessionscombined);
    
end


function [knob_data] = TransformTrialData(files,path)

cd(path);

if ischar(files)                                                        %If only one file was selected...
    files = {files};                                                    %Convert the string to a cell array.
end


knob_data.trial_length = 500;
knob_data.combined_length = 0;
knob_data.num_sessions = length(files);

knob_data.session_length = nan(1, length(files));

knob_data.max_session_trials = 0;

for s = 1:knob_data.num_sessions
    try                                                                     %Try to read in the data file...
        [~,filename,ext] = fileparts(files{s});
        switch ext
            case '.ArdyMotor'
                data = ArdyMotorFileRead(files{s});
                %                 pre_daycodes(f) = data.daycode;
            case '.MotoTrak'
                data = MotoTrakFileRead(files{s});
                data = MotoTrak_to_ArdyMotor(data);
                %                 data.rat = data.subject;
                %                 pre_daycodes(f) = data.start_time;
        end
        %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
    catch err                                                               %If an error occurs...
        warning(['ERROR READING: ' files{1}]);                              %Show which file had a read problem...
        warning(err.message);                                               %Show the actual error message.
    end
%     data = ArdyMotorFileRead(files{s});
    temp = length(data.trial);
    knob_data.combined_length = temp+knob_data.combined_length;
    knob_data.session_length(s) = temp;
    
    if knob_data.session_length(s) > knob_data.max_session_trials;
        knob_data.max_session_trials = knob_data.session_length(s);
    end
end

knob_data.combined_trials = nan(knob_data.combined_length, 500);
knob_data.trial = nan(knob_data.max_session_trials, 500, knob_data.num_sessions);

for s = 1:knob_data.num_sessions
    try                                                                     %Try to read in the data file...
        [~,filename,ext] = fileparts(files{s});
        switch ext
            case '.ArdyMotor'
                data = ArdyMotorFileRead(files{s});
                %                 pre_daycodes(f) = data.daycode;
            case '.MotoTrak'
                data = MotoTrakFileRead(files{s});
                data = MotoTrak_to_ArdyMotor(data);
                %                 data.rat = data.subject;
                %                 pre_daycodes(f) = data.start_time;
        end
        %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
    catch err                                                               %If an error occurs...
        warning(['ERROR READING: ' files{1}]);                              %Show which file had a read problem...
        warning(err.message);                                               %Show the actual error message.
    end
    %     data = ArdyMotorFileRead(files{s});
    for t = 1:knob_data.session_length(s)
        %         if any(data.trial(t).signal < -10)
        %             knob_data.trial(t,:,s) = NaN;
        %             disp('Abnormal negative signal detected.... Removing');
        %         elseif any(data.trial(t).signal > 150)
        %             disp('Abnormal positive signal detected >150 deg... Removing')
        %         else
        knob_data.trial(t, :, s) = data.trial(t).signal(1:500,1);
        %         end
    end
end

for s = 1:knob_data.num_sessions
    for j=1:500
        knob_data.mean_plot(s,j) = nanmean(knob_data.trial(:,j,s));
%         interval(s,j) = nanstd(knob_data.trial(:,j,s))/(sqrt(knob_data.session_length(s)));
        knob_data.interval(s,j) = simple_ci(knob_data.trial(:,j,s));
        knob_data.upper(s,j) = knob_data.mean_plot(s,j) + knob_data.interval(s,j);
        knob_data.lower(s,j) = knob_data.mean_plot(s,j) - knob_data.interval(s,j);
    end
end

function ci = simple_ci(X,alpha)
%Simple confidence interval calculator.

%CI = SIMPLE_CI(X,ALPHA) finds the confidence range for the single column
%dataset X using the significance level ALPHA.  If ALPHA isn't specified,
%the function uses a default value of 0.05.

if (nargin < 2)     %If the use didn't specify an alpha.
    alpha = 0.05; 
end     
if (nargin < 1)
   error('simple_ci requires single column data input.');
end
[h,p,ci] = ttest(X,0,alpha);
ci = nanmean(X,1)-ci(1,:);
