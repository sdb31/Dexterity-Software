function RawWaveformsGUI(vargarin)
[files, path] = uigetfile({'*.MotoTrak;*.ArdyMotor'},...
    'Select Animal(s) Data', ...
    'multiselect','on');
try 
    cd(path);
catch err
    return
end
if ischar(files)                                                      %If only one file was selected...
    files = {files};                                                    %Convert the string to a cell array.
end
for i = 1:length(files)
    D = dir(files{i});
    Bytes(i) = D.bytes;
end
[row, col] = find(Bytes>500);
files = files(col);
knob_data.trial_length = 500;
knob_data.num_sessions = length(files);
knob_data.combined_sessions_length = 0;
knob_data.session_length = nan(1, length(files));
knob_data.max_session_trials = 0;
for s = 1:knob_data.num_sessions
    try                                                                     %Try to read in the data file...
        [~,filename,ext] = fileparts(files{s});
        switch ext
            case '.ArdyMotor'
                data = ArdyMotorFileRead(files{s});
            case '.MotoTrak'
                data = MotoTrakFileRead(files{s});
                data = MotoTrak_to_ArdyMotor(data);
                data.rat = data.subject;
        end
        %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
    catch err                                                               %If an error occurs...
        warning(['ERROR READING: ' files{s}]);                              %Show which file had a read problem...
        warning(err.message);                                               %Show the actual error message.
        continue
    end
%     data = ArdyMotorFileRead(files{s});
    temp = length(data.trial);
    knob_data.combined_length = temp+knob_data.combined_sessions_length;
    knob_data.session_length(s) = temp;
    
    if knob_data.session_length(s) > knob_data.max_session_trials;
        knob_data.max_session_trials = knob_data.session_length(s);
    end
end
knob_data.combined_trials = nan(knob_data.combined_length, 500);
knob_data.trial = nan(knob_data.max_session_trials, 500, knob_data.num_sessions);

%% User Variables
Weeks = inputdlg('How many data sets would you like to analyze?', 'Data Sets', [1 50]);
Weeks = Weeks{:};
Weeks = str2double(Weeks);
Sessions = inputdlg('How many sessions are in each data set?', 'Sessions per Data Set', [1 50]);
Sessions = strsplit(Sessions{:});
Sessions = str2double(Sessions);
Titles = inputdlg('List the titles of each data set followed by a space', 'Data Set Titles', [1 50]);
Titles = strsplit(Titles{:});
Session_Count = zeros(1,length(Sessions)+1);
Waveform_Type = inputdlg('Please write ''Median'' or ''Mean'' (Note: case sensitive)', 'Waveform Analysis', [1 50]);
Waveform_Type = Waveform_Type{:};
%% Initialize Variables
for i = 1:length(Sessions);
    Session_Count(i+1) = sum(Sessions(1:i));
end
for i = 1:length(Sessions);
    Weekly_Trial_Sum(i) = sum(knob_data.session_length((Session_Count(i)+1):(Session_Count(i+1))));
end
Min_Random_Trials = round(0.9*min(Weekly_Trial_Sum));
Max_Distance = nan(2000, length(Sessions));
Peak_Velocity = nan(2000, length(Sessions));
Counter = ones(1,length(Sessions));
Latency_To_Hit = nan(2000, length(Sessions));
for i = 1:length(Weeks);
    knob_data.trial2 = nan(Weekly_Trial_Sum(i), 500, length(Weeks));
    for s = 1:knob_data.num_sessions
        try                                                                     %Try to read in the data file...
            [~,filename,ext] = fileparts(files{s});
            switch ext
                case '.ArdyMotor'
                    data = ArdyMotorFileRead(files{s});
                case '.MotoTrak'
                    data = MotoTrakFileRead(files{s});
                    data = MotoTrak_to_ArdyMotor(data);
                    data.rat = data.subject;
            end
            %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
        catch err                                                               %If an error occurs...
            warning(['ERROR READING: ' files{s}]);                              %Show which file had a read problem...
            warning(err.message);                                               %Show the actual error message.
            continue
        end
        
        %         data = ArdyMotorFileRead(files{s});
        for t = 1:knob_data.session_length(s);
            if size(data.trial(t).signal,1) > size(data.trial(t).signal,2)
                data.trial(t).signal = data.trial(t).signal';
            end
            knob_data.trial(t,:,s) = data.trial(t).signal(1,1:500);
            %             knob_data.trial2(t,i) = data.trial(t).signal;
            knob_data.thresh(t,:,s) = data.trial(t).thresh;
            knob_data.outcome(t,:,s) = data.trial(t).outcome;
        end
    end
end
%% Plots and Analysis
for i = 1:Weeks
    if sum(Sessions(i)) == 1;
        
        figure; clf;
        for m = (Session_Count(i)+1):(Session_Count(i+1));
            try                                                                     %Try to read in the data file...
                [~,filename,ext] = fileparts(files{m});
                switch ext
                    case '.ArdyMotor'
                        data = ArdyMotorFileRead(files{m});
                    case '.MotoTrak'
                        data = MotoTrakFileRead(files{m});
                        data = MotoTrak_to_ArdyMotor(data);
                        data.rat = data.subject;
                end
                %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
            catch err                                                               %If an error occurs...
                warning(['ERROR READING: ' files{m}]);                              %Show which file had a read problem...
                warning(err.message);                                               %Show the actual error message.
                continue
            end
%             data = ArdyMotorFileRead(files{m});
            for t = 1:knob_data.session_length(m)
                TempMatrix(Counter(i),:) = data.trial(t).signal(1:500,1)';
                Max_Distance(Counter(i),i) = max(knob_data.trial(t,:,m));
                a = find((data.trial(t).sample_times < 1000*data.trial(t).hitwin));
                if exist('data.threshtype','var') == 0;
                    data.threshtype = [];
                end
                if strcmpi(data.threshtype,'degrees (bidirectional)')               %If the threshold type is bidirectional knob-turning...
                    signal = abs(data.trial(t).signal(a) - ....
                        data.trial(t).signal(1));                                   %Subtract the starting degrees value from the trial signal.
                elseif strcmpi(data.threshtype,'# of spins')                        %If the threshold type is the number of spins...
                    temp = diff(data.trial(t).signal);                              %Find the velocity profile for this trial.
                    temp = boxsmooth(temp,10);                                      %Boxsmooth the wheel velocity with a 100 ms smooth.
                    [pks,i] = PeakFinder(temp,10);                                  %Find all peaks in the trial signal at least 100 ms apart.
                    i(pks < 1) = [];                                                %Kick out all peaks that are less than 1 degree/sample.
                    i = intersect(a,i+1)-1;                                         %Find all peaks that are in the hit window.
                    signal = length(i);                                             %Set the trial signal to the number of wheel spins.
                elseif strcmpi(data.threshtype, 'degrees (total)')
                    signal = data.trial(t).signal(a);               
                else
                    signal = data.trial(t).signal(a) - data.trial(t).signal(1);    %Grab the raw signal for the trial.
                end
                smooth_trial_signal = boxsmooth(signal);
                Peak_Velocity(Counter(i),i) = max(boxsmooth(diff(smooth_trial_signal)))*100;
                if (knob_data.outcome(t,:,m) == 72)                                   %If it was a hit
                    hit_time = find(knob_data.trial(t,:,m) >= ...                    %Calculate the hit time
                        knob_data.thresh(t,:,m),1);
                    Latency_To_Hit(Counter(i),i) = hit_time;                       %hit time is then latency to hit
                else
                    Latency_To_Hit(Counter(i),i) = NaN;                            %If trial resulted in a miss, then set latency to hit to NaN
                end
                Counter(i) = Counter(i) + 1;
            end
            for j = 1:500;
                Mean_Plot(m,j) = nanmean(knob_data.trial(:,j,m));
                Median_Plot(m,j) = nanmedian(knob_data.trial(:,j,m));
                SD_Plot(m,j) = nanstd(knob_data.trial(:,j,m));
            end
        end
        TempMatrix = datasample(TempMatrix,Min_Random_Trials,'Replace', false);
        hold on;
        for t = 1:Min_Random_Trials;
            patchline(1:500, TempMatrix(t,:), 'edgecolor', 'b', 'linewidth', .5, 'edgealpha', 0.075);
        end
        hold off;
        Overall_Mean(i,:) = Mean_Plot((Session_Count(i)+1):(Session_Count(i+1)),:);
        Overall_Median(i,:) = Median_Plot((Session_Count(i)+1):(Session_Count(i+1)),:);
        Overall_SD(i,:) = SD_Plot((Session_Count(i)+1):(Session_Count(i+1)),:);
        title(Titles(i), 'Fontsize', 10, 'Fontweight', 'normal');
        hold on;
        boxplot(Max_Distance(:,i), 'positions', 40, 'widths', 50, 'outliersize', 6, 'colors', 'b', 'symbol', 'b.');
        hold off;
        hold on;
        boxplot(Latency_To_Hit(:,i), 'orientation', 'horizontal', 'widths', 30, 'positions', -20, 'outliersize', 6, 'colors', 'b', 'symbol', 'b.');
        hold off;
        YTickLabels = -40:5:round(max(Max_Distance(:,i))/5)*5;
        %                 set(gca, 'TickDir', 'out', 'YLim', [-40 round(1.1*max(Max_Distance(:,i)))], 'YTick', -40:40:round(1.1*max(Max_Distance(:,i))), 'XLim', [0 350], 'XTick', 0:50:350, 'XTickLabels', {'0', '50', '100', '150', '200', '250', '300', '350'},...
        %                     'YTickLabels', YTickLabels);
        set(gca, 'TickDir', 'out', 'YLim', [-40 round(1.1*max(Max_Distance(:,i)))],...
            'YTick',-40:5:round(max(Max_Distance(:,i))/5)*5,'YTickLabels', YTickLabels,...
            'XLim', [0 350], 'XTick', 0:50:350, 'XTickLabels', {'0', '50', '100', '150', '200', '250', '300', '350'});
        %         ylabel('Angle (degrees)', 'Fontsize', 10);
        box off;
    elseif sum(Sessions(i)) > 1;
        
        figure; clf;
        TempMatrix = [];
        for m = (Session_Count(i)+1):(Session_Count(i+1));
            try                                                                     %Try to read in the data file...
                [~,filename,ext] = fileparts(files{m});
                switch ext
                    case '.ArdyMotor'
                        data = ArdyMotorFileRead(files{m});
                    case '.MotoTrak'
                        data = MotoTrakFileRead(files{m});
                        data = MotoTrak_to_ArdyMotor(data);
                        data.rat = data.subject;
                end
                %         temp = ArdyMotorFileRead(files{f});                                 %Read in the data from each file.
            catch err                                                               %If an error occurs...
                warning(['ERROR READING: ' files{m}]);                              %Show which file had a read problem...
                warning(err.message);                                               %Show the actual error message.
                continue
            end
            %             data = ArdyMotorFileRead(files{m});
            for t = 1:knob_data.session_length(m);
                TempMatrix(Counter(i),:) = data.trial(t).signal(1:500,1)';
                Max_Distance(Counter(i),i) = max(knob_data.trial(t,:,m));
                a = find((data.trial(t).sample_times < 1000*data.trial(t).hitwin));
                if exist('data.threshtype','var') == 0;
                    data.threshtype = [];
                end
                if strcmpi(data.threshtype,'degrees (bidirectional)')               %If the threshold type is bidirectional knob-turning...
                    signal = abs(data.trial(t).signal(a) - ....
                        data.trial(t).signal(1));                                   %Subtract the starting degrees value from the trial signal.
                elseif strcmpi(data.threshtype,'# of spins')                        %If the threshold type is the number of spins...
                    temp = diff(data.trial(t).signal);                              %Find the velocity profile for this trial.
                    temp = boxsmooth(temp,10);                                      %Boxsmooth the wheel velocity with a 100 ms smooth.
                    [pks,i] = PeakFinder(temp,10);                                  %Find all peaks in the trial signal at least 100 ms apart.
                    i(pks < 1) = [];                                                %Kick out all peaks that are less than 1 degree/sample.
                    i = intersect(a,i+1)-1;                                         %Find all peaks that are in the hit window.
                    signal = length(i);                                             %Set the trial signal to the number of wheel spins.
                elseif strcmpi(data.threshtype, 'degrees (total)')
                    signal = data.trial(t).signal(a);
                else
                    signal = data.trial(t).signal(a) - data.trial(t).signal(1);    %Grab the raw signal for the trial.
                end
                smooth_trial_signal = boxsmooth(signal);
                Peak_Velocity(Counter(i),i) = max(boxsmooth(diff(smooth_trial_signal)))*100;
                if (knob_data.outcome(t,:,m) == 72)                                   %If it was a hit
                    hit_time = find(knob_data.trial(t,:,m) >= ...                    %Calculate the hit time
                        knob_data.thresh(t,:,m),1);
                    Latency_To_Hit(Counter(i),i) = hit_time;                       %hit time is then latency to hit
                else
                    Latency_To_Hit(Counter(i),i) = NaN;                            %If trial resulted in a miss, then set latency to hit to NaN
                end
                Counter(i) = Counter(i) + 1;
            end
            for j = 1:500;
                Mean_Plot(m,j) = nanmean(knob_data.trial(:,j,m));
                Median_Plot(m,j) = nanmedian(knob_data.trial(:,j,m));
                SD_Plot(m,j) = nanstd(knob_data.trial(:,j,m));
            end
        end
        TempMatrix = datasample(TempMatrix,Min_Random_Trials,'Replace', false);
        hold on;
        for t = 1:Min_Random_Trials;
            patchline(1:500, TempMatrix(t,:), 'edgecolor', 'b', 'linewidth', .5, 'edgealpha', 0.05);
        end
        hold off;
        Overall_Mean(i,:) = nanmean(Mean_Plot((Session_Count(i)+1):(Session_Count(i+1)),:));
        Overall_Median(i,:) = nanmedian(Median_Plot((Session_Count(i)+1):(Session_Count(i+1)),:));
        Overall_SD(i,:) = nanmean(SD_Plot((Session_Count(i)+1):(Session_Count(i+1)),:));
        title(Titles(i), 'Fontsize', 10, 'Fontweight', 'normal');
        hold on;
        boxplot(Max_Distance(:,i), 'positions', 40, 'widths', 50, 'outliersize', 6, 'colors', 'b', 'symbol', 'b.');
        hold off;
        hold on;
        boxplot(Latency_To_Hit(:,i), 'orientation', 'horizontal', 'widths', 30, 'positions', -20, 'outliersize', 6, 'colors', 'b', 'symbol', 'b.');
        hold off;
        YTickLabels = -40:5:round(max(Max_Distance(:,i))/5)*5;
        %                 set(gca, 'TickDir', 'out', 'YLim', [-40 round(1.1*max(Max_Distance(:,i)))], 'YTick', -40:40:round(1.1*max(Max_Distance(:,i))), 'XLim', [0 350], 'XTick', 0:50:350, 'XTickLabels', {'0', '50', '100', '150', '200', '250', '300', '350'},...
        %                     'YTickLabels', YTickLabels);
        set(gca, 'TickDir', 'out', 'YLim', [-40 round(1.1*max(Max_Distance(:,i)))],...
            'YTick',-40:5:round(max(Max_Distance(:,i))/5)*5,'YTickLabels', YTickLabels,...
            'XLim', [0 350], 'XTick', 0:50:350, 'XTickLabels', {'0', '50', '100', '150', '200', '250', '300', '350'});
        %         ylabel('Angle (degrees)', 'Fontsize', 10);
        box off;
    end
end


%% Average Waveforms Plots
ThreeDee = questdlg('Do you want to plot Mean or Median Waveform in 3D?',...
    '3D Option',...
    'Yes','No','No');
Shading_Option = questdlg('Do you want to include error shading?',...
    'Error Shading',...
    'SEM','None','None');
FF = Weeks + 2;
figure; clf;
switch Waveform_Type
    case 'Median'
        for d = 1:size(Overall_Median,1)
            if strcmpi(ThreeDee,'Yes') == 1;
                hold on;
                if strcmpi(Shading_Option,'SEM') == 1;
                    SEM_Upper = Overall_Median(d,:) + Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    SEM_Lower = Overall_Median(d,:) - Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    fill3([1:500,500:-1:1], ...
                        0.5*d*ones(1,1000),...
                        [SEM_Upper, ...
                        fliplr(SEM_Lower)], ...
                        [0 0.5 0]);
                end
                h(d) = plot3( 1:500, 0.5*d*ones(1,500),Overall_Median(d,:),'Linewidth',2)
                hold off;
                grid on; xlabel('Time (hs)'); zlabel('Dependent Variable');
                title('Median Waveforms', 'Fontweight', 'Normal', 'Fontsize', 10);
                legend(h, Titles);
            else
                hold on;
                if strcmpi(Shading_Option,'SEM') == 1;
                    SEM_Upper = Overall_Median(d,:) + Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    SEM_Lower = Overall_Median(d,:) - Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    fill([1:500,500:-1:1], ...
                        [SEM_Upper, ...
                        fliplr(SEM_Lower)], ...
                        [0 0.5 0]);
                end
                h(d) = plot(1:500, Overall_Median(d,:), 'Linewidth', 2);
                hold off;
                title('Median Waveforms', 'Fontweight', 'Normal', 'Fontsize', 10);
                ylabel('Angle (degrees)', 'Fontsize', 10);
                legend(h, Titles);
            end
        end
    case 'Mean'
        for d = 1:size(Overall_Mean,1)
            if strcmpi(ThreeDee,'Yes') == 1;
                hold on;
                if strcmpi(Shading_Option,'SEM') == 1;
                    SEM_Upper = Overall_Mean(d,:) + Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    SEM_Lower = Overall_Mean(d,:) - Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    fill3([1:500,500:-1:1], ...
                        0.5*d*ones(1,1000),...
                        [SEM_Upper, ...
                        fliplr(SEM_Lower)], ...
                        [0 0.5 0]);
                end
                h(d) = plot3( 1:500, 0.5*d*ones(1,500),Overall_Mean(d,:),'Linewidth',2);
                hold off;
                grid on; xlabel('Time (hs)'); zlabel('Dependent Variable');
                title('Mean Waveforms', 'Fontweight', 'Normal', 'Fontsize', 10);
                legend(h,Titles);
            else
                hold on;
                if strcmpi(Shading_Option,'SEM') == 1;
                    %                     legend(Titles);
                    SEM_Upper = Overall_Mean(d,:) + Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    SEM_Lower = Overall_Mean(d,:) - Overall_SD(d,:)/sqrt(Min_Random_Trials);
                    fill([1:500,500:-1:1], ...
                        [SEM_Upper, ...
                        fliplr(SEM_Lower)], ...
                        [0 0.5 0]);
                end
                h(d) = plot(1:500, Overall_Mean(d,:), 'Linewidth', 2);
                hold off
                title('Mean Waveforms', 'Fontweight', 'Normal', 'Fontsize', 10);
                ylabel('Angle (degrees)', 'Fontsize', 10);
                legend(h,Titles);
            end
        end

end
% legend(Titles);
box off;
set(gca, 'TickDir', 'out');

figure; clf;
boxplot(Peak_Velocity, 'symbol', 'k.', 'outliersize', 3, 'colors', 'k');
set(gca, 'TickDir', 'out', 'XTickLabels', Titles);
box off;
title('Peak Velocity', 'Fontweight', 'normal', 'Fontsize', 10);
ylabel('Velocity (deg/s)', 'Fontsize', 10);

figure; clf;
boxplot(Max_Distance, 'symbol', 'k.', 'outliersize', 3, 'colors', 'k');
set(gca, 'TickDir', 'out', 'XTickLabels', Titles);
box off;
title('Signal Peak', 'Fontweight', 'normal', 'Fontsize', 10);
ylabel('Dependent Variable', 'Fontsize', 10);

figure; clf;
boxplot(Latency_To_Hit, 'symbol', 'k.', 'outliersize', 3, 'colors', 'k');
set(gca, 'TickDir', 'out', 'XTickLabels', Titles);
box off;
title('Latency to HIT', 'Fontweight', 'normal', 'Fontsize', 10);
ylabel('Latency (hs)', 'Fontsize', 10);
end

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
end