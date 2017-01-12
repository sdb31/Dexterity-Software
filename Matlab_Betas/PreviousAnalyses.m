function PreviousAnalyses(varargin)
datapath = 'C:\KnobAnalysis\ConfigFiles\';                                         %Set the primary local data path for saving data files.
if ~exist(datapath,'dir')                                           %If the primary local data path doesn't already exist...
    mkdir(datapath);                                                %Make the primary local data path.
end
path = 'C:\AnalyzeGroup';    

Info = dir(path);
ExperimentNames = {Info.name};
ExperimentNames = ExperimentNames(3:end);
if isempty(ExperimentNames) ==1;
    ExperimentNames = {'No Analyses'};
end

Info = dir(datapath);
AnalysisNames = {Info.name}; AnalysisNames = AnalysisNames(3:end);
if isempty(AnalysisNames) ==1;
    AnalysisNames = {'No Analyses'};
end

set(0,'units','centimeters');                                               %Set the system units to centimeters.
pos = get(0,'Screensize');                                                  %Grab the screensize.
h = 8;                                                                      %Set the height of the figure.
w = 15;                                                                     %Set the width of the figure.
handles.fig = figure('numbertitle','off','name','Dexterity: Previous Analyses',...
    'units','centimeters','Position',[pos(3)/2-w/2, pos(4)/2-h/2, w, h],...
    'menubar','none','resize','off');
handles.Non_Annotated = uipanel('Parent',handles.fig,'Title','Non-Annotated Analyses',...
    'Fontsize',12,'Fontweight','bold',...
    'Position',[.025 .05 .45 .9]);
handles.Non_Annotated_List = uicontrol('Parent',handles.Non_Annotated,'style', 'listbox','HorizontalAlignment', 'left',...
    'units','normalized','position',[.05 .05 .9 .9],'string',AnalysisNames,...
    'Fontsize',10);
handles.Annotated = uipanel('Parent',handles.fig,'Title','Annotated Analyses',...
    'Fontsize',12,'Fontweight','bold',...
    'Position',[.525 .05 .45 .9]);
handles.Annotated_List = uicontrol('Parent',handles.Annotated,'style', 'listbox','HorizontalAlignment', 'left',...
    'units','normalized','position',[.05 .05 .9 .9],'string',ExperimentNames,...
    'Fontsize',10);

set(handles.Non_Annotated_List, 'callback', {@ExistingAnalysis_Jump});
set(handles.Annotated_List, 'callback', {@ExistingAnnotatedAnalysis_Jump});


end

function ExistingAnalysis_Jump(hObject,~)
Existing_Analysis(hObject)
end

function ExistingAnnotatedAnalysis_Jump(hObject,~)
ExistingAnnotatedAnalysis(hObject)
end