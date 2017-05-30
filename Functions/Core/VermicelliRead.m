function [data] = VermicelliRead(files)
for i = 1:length(files)                                                                                                 % Step through each Vermicelli Master File
    [~,sheets] = xlsfinfo(char(files(i)));                                                                         % Read in the info for each sheet within each Master FIle
    No_Sheets = numel(sheets);                                                                                          % Calculate the number of sheets
    if No_Sheets > 1;                                                                                                   % If the number of sheets is greater than 1 AKA the number of animals is greater than 1
        for s = 1:No_Sheets;                                                                                            % Step through each sheet
            [NUM,TXT,~] = xlsread(char(files(i)),s);                                                                  % Read in number, text, and all raw data from each sheet
            [~,~,~] = xlsfinfo(char(files(i)));
            data(i).animal(s).name = TXT(2,2);                                                                          % Read in the animal name from the second row, second column of the TXT matrix
            for q = 1:(size(NUM,2)/3);                                                                                  % Step through each block of columns per date
                data(i).animal(s).date{q} = TXT{3,3*q-1};                                                               % Read in the date                
                data(i).animal(s).video{q} = NUM(1:end,3*q-2)';                                                         % Read in the number of videos
                data(i).animal(s).targeted{q} = NUM(1:end,3*q-1)';                                                      % Read in the number of targeted manipulations
                data(i).animal(s).nontargeted{q} = NUM(1:end,3*q)';                                                     % Read in the number of non-targeted manipulations
            end
            for l = 1:(size(data(i).animal(s).targeted,2));                                                             % Calculate average stats for each day
                data(i).animal(s).mean_targeted(l) = nanmean(data(i).animal(s).targeted{l});                            % Calculate mean number of targeted manipulations
                data(i).animal(s).mean_nontargeted(l) = nanmean(data(i).animal(s).nontargeted{l});                      % Calculate mean number of non-targeted manipulations
            end
            data(i).animal(s).ratio = log2(data(i).animal(s).mean_targeted./data(i).animal(s).mean_nontargeted);        % Calculate the ratio of targeted to non-targeted manipulations
            clear NUM TXT RAW                                                                                           % Clear the variables before stepping through the next sheet
        end
    else                                                                                                                % If the number of sheets isn't greater than 1 AKA the number of animals equals 1...
        [NUM,TXT,~] = xlsread(char(files(i)));                                                                        % ...repeat the steps seen above for 1 animal
        data(i).animal.name = TXT(2,2);
        for q = 1:(size(NUM,2)/3);
            data(i).animal.date{q} = TXT{3,3*q-1};
            data(i).animal.video{q} = NUM(1:end,3*q-2)';
            data(i).animal.targeted{q} = NUM(1:end,3*q-1)';
            data(i).animal.nontargeted{q} = NUM(1:end,3*q)';
        end
        for l = 1:(size(data(i).animal(s).targeted,2));
            data(i).animal.mean_targeted(l) = nanmean(data(i).animal(s).targeted{l});
            data(i).animal.mean_nontargeted(l) = nanmean(data(i).animal(s).nontargeted{l});
        end
        data(i).animal.ratio = log2(data(i).animal(s).mean_targeted./data(i).animal(s).mean_nontargeted);
        clear NUM TXT RAW
    end
    data(i).filename = files(i);                                                                                        % Save the Filename
    data(i).device = 'Vermicelli';                                                                                      % Save the device name. Dexterity requires this for analysis
    clear status sheets No_Sheets
end