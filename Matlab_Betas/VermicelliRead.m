function [data] = VermicelliRead(files);
% datapath = 'C:\';
% datapath = uigetdir(datapath, 'Where is your Vermicelli data located?');       %Ask the user where their data is located.
% if datapath(1) == 0                                                         %If the user pressed "cancel"...
%     return                                                                  %Skip execution of the rest of the function.
% end
% 
% files = file_miner(datapath,{'*.xls','*.xlsx'});  
% pause(0.01);                                                                %Pause for 10 milliseconds.
% if isempty(files)                                                           %If no files were found...
%     errordlg('No Vermicelli files were found!');   %Show an error dialog box.
% end
% % 
% for i = 1:length(files)
%     [status,sheets] = xlsfinfo(char(files(i)));
%     No_Sheets = numel(sheets);
%     if No_Sheets > 1;
%         for s = 1:No_Sheets;
%            T = readtable(char(files(i)),'Sheet', s);         
%            data(i).animal(s).name = T{1,2};
%            for q = 1:(size(T,2)/3);
%                data(i).animal(s).video{q} = T{4:end,3*q-2}';  
%                data(i).animal(s).targeted{q} = T{4:end,3*q-1}';  
%                data(i).animal(s).nontargeted{q} = T{4:end,3*q}';  
%            end
%            clear T
%         end
%     else
%         T = readtable(char(files(i)));
%         data(i).animal(1).name = T{1,2};
%         clear T
%     end
%     data(i).filename = files(i);
%     clear status sheets No_Sheets
% end

for i = 1:length(files)
    [status,sheets] = xlsfinfo(char(files(i)));
    No_Sheets = numel(sheets);
    if No_Sheets > 1;
        for s = 1:No_Sheets;
            [NUM,TXT,RAW] = xlsread(char(files(i)),s);
            [status,sheets,xlFormat] = xlsfinfo(char(files(i)));
            data(i).animal(s).name = TXT(2,2);
            for q = 1:(size(NUM,2)/3);
                data(i).animal(s).date{q} = TXT{3,3*q-1};
                data(i).animal(s).video{q} = NUM(1:end,3*q-2)';
                data(i).animal(s).targeted{q} = NUM(1:end,3*q-1)';
                data(i).animal(s).nontargeted{q} = NUM(1:end,3*q)';
            end
            for l = 1:(size(data(i).animal(s).targeted,2));
                data(i).animal(s).mean_targeted(l) = nanmean(data(i).animal(s).targeted{l});
                data(i).animal(s).mean_nontargeted(l) = nanmean(data(i).animal(s).nontargeted{l});                
            end
            data(i).animal(s).ratio = log2(data(i).animal(s).mean_targeted./data(i).animal(s).mean_nontargeted);
            clear NUM TXT RAW
        end
    else
        [NUM,TXT,RAW] = xlsread(char(files(i)));
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
    data(i).filename = files(i);
    data(i).device = 'Vermicelli';
    clear status sheets No_Sheets
end