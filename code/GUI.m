function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 16-Oct-2015 10:14:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.SelectFileButton,'TooltipString', sprintf(help_msg('raw_data')))
set(handles.ConsideredInputListbox,'TooltipString', sprintf(help_msg('select_io')))
set(handles.ConsideredOutputListbox,'TooltipString', sprintf(help_msg('select_io')))
set(handles.TimeVariableCheckbox,'TooltipString', sprintf(help_msg('include_time')))
set(handles.kfoldTextbox,'TooltipString', sprintf(help_msg('folds')))
set(handles.PermutationTextbox,'TooltipString', sprintf(help_msg('permutations')))
%set(handles.open_final_data,'TooltipString', sprintf(help_msg('final_data')))
%set(handles.save_final_data,'TooltipString', sprintf(help_msg('final_data')))
set(handles.showdata,'Enable','off')
set(handles.save_analysis,'Enable','off')
set(handles.popup_xy,'Enable','off');
set(handles.show_plots,'Enable','off');
set(handles.display_analysis,'Enable','off');
end


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Select your file.
function SelectFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname,~] = uigetfile({'*.csv';'*.xls';'*.xlsx';'*.mat'},'Select a valid input file','../raw_data_csv');

if ischar(filename) && ischar(pathname)
    
    error = 0;
    
    if strcmp(filename(end-2:end),'mat')
        %% MAT-file import
        load([pathname '/' filename]);
        handles.in.name = filename(1:end-4);
    elseif strcmp(filename(end-2:end),'csv') || strcmp(filename(end-2:end),'xls') || strcmp(filename(end-3:end),'xlsx')
        %% XLS/CSV/XLSX-file import
        set(handles.ErrorDisplay,'String','')
        if strcmp(filename(end-3:end),'xlsx') || strcmp(filename(end-2:end),'xls')
            [data,header,~] = xlsread([pathname '/' filename]);
        elseif strcmp(filename(end-2:end),'csv')
            [header,data] = csvreadh([pathname '/' filename]);
        end
        if strcmp(filename(end-2:end),'csv') || strcmp(filename(end-2:end),'xls')
            handles.in.name = filename(1:end-4);
        elseif strcmp(filename(end-3:end),'xlsx')
            handles.in.name = filename(1:end-5);
        end
        if length(header)<1
            for i = 1:size(data,2)
                header{1,i} = ['Var' num2str(i)];
            end
        end
    else
        error = 1;
        ResetScript
        set(handles.ErrorDisplay,'String','Wrong input file extension.')
        set(handles.ErrorDisplay,'ForegroundColor','r')
    end

    if error == 0

        set(handles.ErrorDisplay,'String','')

        handles.in.data = data;
        handles.n = size(handles.in.data,1);
        handles.in.filename = filename;
        if ~exist('header','var')%if no variable called header, we create it
            for i = 1:size(handles.in.data,2)
                header{1,i} = ['Var' num2str(i)];
                handles.in.headers = header;
            end
        end

        [m_in,n_in] = size(header);

        %Check whether we have a column or line vector (line vector
        %required)
        if m_in < n_in
            handles.in.headers = header';
        else
            handles.in.headers = header;
        end

        set(handles.ConsideredInputListbox,'Value',[1]);
        set(handles.ConsideredOutputListbox,'Value',[1]);
        set(handles.ConsideredInputListbox,'Enable','on')
        set(handles.ConsideredInputListbox,'String',handles.in.headers)
        set(handles.ConsideredInputListbox,'Max',max(m_in,n_in))
        set(handles.ConsideredOutputListbox,'Enable','on')
        set(handles.ConsideredOutputListbox,'String',handles.in.headers)
        set(handles.ConsideredOutputListbox,'Max',1) %Only one single output can be selected at a time
        set(handles.description,'String',' ')
        set(handles.TimeVariableCheckbox,'Enable','on')
        set(handles.kfoldTextbox,'Enable','on')
        set(handles.PermutationTextbox,'Enable','on')
        set(handles.ngrid,'Enable','on')
        set(handles.RunButton,'Enable','on')
        set(handles.showdata,'Enable','on');
        set(handles.popup_xy,'Enable','off');
        set(handles.show_plots,'Enable','off');
        set(handles.display_analysis,'Enable','off');
        set(handles.save_analysis,'Enable','off');
        cla(handles.MainPlot);
        cla(handles.SecondPlot);
        cla(handles.ThirdPlot);
    end
    contents = cellstr(get(handles.ConsideredOutputListbox,'String'));
    index_sel = get(handles.ConsideredOutputListbox,'Value');
    handles.in.considered_output=contents(index_sel');
    contents = cellstr(get(handles.ConsideredInputListbox,'String'));
    index_sel = get(handles.ConsideredInputListbox,'Value');
    handles.in.considered_inputs=contents(index_sel');
guidata(hObject,handles);
end
end

% --- Executes on selection change in ConsideredInputListbox.
function ConsideredInputListbox_Callback(hObject, eventdata, handles)
% hObject    handle to ConsideredInputListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ConsideredInputListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ConsideredInputListbox

contents = cellstr(get(hObject,'String'));
index_sel = get(hObject,'Value');
handles.in.considered_inputs=contents(index_sel');
%handles.in.considered_input_index = index_sel;

% Clear variables that depend on in the selected inputs:
if isfield(handles.in,'hyp')
    handles.in = rmfield(handles.in,'hyp');
end
if isfield(handles.in,'covfunction')
    handles.in = rmfield(handles.in,'covfunction');
end
if isfield(handles.in,'plot_xy')
    handles.in = rmfield(handles.in,'plot_xy');
end
if isfield(handles.in,'idx_xy')
    handles.in = rmfield(handles.in,'idx_xy');
end
set(handles.save_analysis,'Enable','off');

guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function ConsideredInputListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ConsideredInputListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in ConsideredOutputListbox.
function ConsideredOutputListbox_Callback(hObject, eventdata, handles)
% hObject    handle to ConsideredOutputListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ConsideredOutputListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ConsideredOutputListbox

contents = cellstr(get(hObject,'String'));
index_sel = get(hObject,'Value');
handles.in.considered_output=contents(index_sel');
%handles.in.considered_output_index = index_sel;

% Clear variables that depend on in the selected inputs:
if isfield(handles.in,'hyp')
    handles.in = rmfield(handles.in,'hyp');
end
if isfield(handles.in,'covfunction')
    handles.in = rmfield(handles.in,'covfunction');
end
set(handles.save_analysis,'Enable','off');


guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function ConsideredOutputListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ConsideredOutputListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in TimeVariableCheckbox.
function TimeVariableCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to TimeVariableCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TimeVariableCheckbox

handles.in.consider_temporality = get(hObject,'Value');
set(handles.save_analysis,'Enable','off');

guidata(hObject,handles);
end



function kfoldTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to kfoldTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kfoldTextbox as text
%        str2double(get(hObject,'String')) returns contents of kfoldTextbox as a double

handles.in.kfolds = str2double(get(hObject,'String'));

set(handles.save_analysis,'Enable','off');

guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function kfoldTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kfoldTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function PermutationTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to PermutationTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PermutationTextbox as text
%        str2double(get(hObject,'String')) returns contents of PermutationTextbox as a double

handles.in.perm = str2double(get(hObject,'String'));

set(handles.save_analysis,'Enable','off');

guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function PermutationTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PermutationTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.MainPlot)
cla(handles.SecondPlot)
cla(handles.ThirdPlot);
set( gcf, 'toolbar', 'none' )

% Get all values from the different form components:
handles.in.kfolds = str2double(get(handles.kfoldTextbox,'String'));
handles.in.perm = str2double(get(handles.PermutationTextbox,'String'));
handles.in.Ngrid = str2double(get(handles.ngrid,'String'));
set(handles.ErrorDisplay,'String','');
set(handles.ErrorDisplay,'ForegroundColor','k');

[handles.in,errors,warnings] = inputs_sanity_check(handles.in);

if ~isempty(errors)
    set(handles.ErrorDisplay,'String',errors{1});
    set(handles.ErrorDisplay,'ForegroundColor','r');
else
    set(handles.ErrorDisplay,'String','Analysis started, check the Matlab Console');
    set(handles.ErrorDisplay,'ForegroundColor','b');
    set(hObject,'Enable','off');
    set(handles.figure1, 'pointer', 'watch');
    drawnow;

    handles.results = main_model(handles.in);
    
    %set(hObject,'Enable','on')
    set(handles.ErrorDisplay,'String','Analysis terminated');
    set(handles.ErrorDisplay,'ForegroundColor','b');   
    set(handles.figure1, 'pointer', 'arrow')
    set(handles.save_analysis,'Enable','on')
    set(handles.display_analysis,'Enable','on')
    set(handles.show_plots,'Enable','on')
    set(hObject,'Enable','on');
    
    Ninputs = length(handles.in.considered_inputs);
    if Ninputs > 2
        comb = nchoosek(1:Ninputs,2);
        liste = {'Most relevant'};
        for i = 1:size(comb,1)
            liste{i+1} = [handles.in.considered_inputs{comb(i,1)} ' - ' handles.in.considered_inputs{comb(i,2)}];
        end
        set(handles.popup_xy,'String',liste)
        set(handles.popup_xy,'Enable','on')
    else
        set(handles.popup_xy,'String',{'Most relevant'})
        set(handles.popup_xy,'Enable','off')
    end
    drawnow;
    
    % Plotting:
    cla(handles.MainPlot);
    cla(handles.SecondPlot);
    cla(handles.ThirdPlot);
    axes(handles.MainPlot);
    plot_regression(handles.in,handles.results);
    axes(handles.SecondPlot);
    plot_prediction(handles.in,handles.results);
    axes(handles.ThirdPlot);
    plot_gaussian(handles.in,handles.results);
    
    %% HTML report

    %Report construction
    % the buildreport does not seem to work in Matlab 2014. Desactivated
    % until it can be debugged.
    %report(buildReport(handles.in,handles.results));

    %Suppression of the image createn for the report
    
    %delete([handles.in.name '_fig1.png']);
    %delete([handles.in.name '_fig2.png']);
    %delete([handles.in.name '_fig3.png']);
    
    result_analysis(handles.in,handles.results)
end
guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function MainPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MainPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate MainPlot
end

% --- Executes on button press in load_analysis.
function load_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to load_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile({'*.mat'},'Select a valid simulation result file','../saved_simulations');
fullpath = [path file];
load(fullpath);

if ~exist('in') || ~exist('results')
    set(handles.ErrorDisplay,'String','The .mat file could not be loaded because it does not contain the "in" and "results" variables. It might not be a proper GPExp data file');
    set(handles.ErrorDisplay,'ForegroundColor','r');
else
    handles.in = in;
    handles.results = results;
    set(handles.ConsideredInputListbox,'Enable','on')
    set(handles.ConsideredInputListbox,'String',handles.in.headers)
    set(handles.ConsideredInputListbox,'Max',length(handles.in.headers))
    set(handles.ConsideredOutputListbox,'Enable','on')
    set(handles.ConsideredOutputListbox,'String',handles.in.headers)
    set(handles.ConsideredOutputListbox,'Max',1) %Only one single output can be selected at a time
    set(handles.TimeVariableCheckbox,'Value',handles.in.consider_temporality)   
    set(handles.kfoldTextbox,'String',num2str(handles.in.kfolds))
    set(handles.PermutationTextbox,'String',num2str(handles.in.perm))
    set(handles.description,'String',handles.in.description)
    set(handles.ngrid,'String',num2str(handles.in.Ngrid))
    set(handles.figure1, 'pointer', 'arrow')
    set(handles.ErrorDisplay,'String',['Loaded previous simulation: ' handles.in.name]);
    set(handles.ErrorDisplay,'ForegroundColor','b');  
    
    set(handles.PermutationTextbox,'Enable','on')    
    set(handles.ngrid,'Enable','on')  
    set(handles.RunButton,'Enable','on')
    set(handles.showdata,'Enable','on');
    set(handles.kfoldTextbox,'Enable','on')
    set(handles.showdata,'Enable','on');
    set(handles.show_plots,'Enable','on');
    set(handles.save_analysis,'Enable','off')
    set(handles.display_analysis,'Enable','on')
    
    Ninputs = length(handles.in.considered_inputs);
    idx = zeros(Ninputs,1);
    for i = 1:Ninputs
        idx(i) = find(strcmp([handles.in.headers], handles.in.considered_inputs{i}));
    end
    idx_out = find(strcmp([handles.in.headers], handles.in.considered_output{1}));
    set(handles.ConsideredInputListbox,'Value',idx);
    set(handles.ConsideredOutputListbox,'Value',idx_out);
    
    if Ninputs > 2
        comb = nchoosek(1:Ninputs,2);
        liste = {'Most relevant'};
        for i = 1:size(comb,1)
            liste{i+1} = [handles.in.considered_inputs{comb(i,1)} ' - ' handles.in.considered_inputs{comb(i,2)}];
        end
        set(handles.popup_xy,'String',liste)
        set(handles.popup_xy,'Enable','on')
    else
        set(handles.popup_xy,'String',{'Most relevant'})
        set(handles.popup_xy,'Enable','off')
    end
    
    % Plotting:
    cla(handles.MainPlot,'reset')
    axes(handles.MainPlot);
    plot_regression(handles.in,handles.results);
    cla(handles.SecondPlot,'reset')
    axes(handles.SecondPlot);
    plot_prediction(handles.in,handles.results);
    cla(handles.ThirdPlot,'reset')
    axes(handles.ThirdPlot);
    plot_gaussian(handles.in,handles.results);
end

guidata(hObject,handles);
end

% --- Executes on button press in save_analysis.
function save_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to save_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'in') || ~isfield(handles,'results')
    set(handles.ErrorDisplay,'String','Results do no seem to have been computed');
    set(handles.ErrorDisplay,'ForegroundColor','r');
    set(hObject,'Enable','off')
else
    in = handles.in;
    results = handles.results;
    [file,path] = uiputfile('*.mat','Save simulation','../saved_simulations');
    fullpath = [path file];
    save(fullpath,'in','results');
end
    
end

% --- Executes on button press in show_plots.
function show_plots_Callback(hObject, eventdata, handles)
% hObject    handle to show_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Plotting:
figure(1)
plot_regression(handles.in,handles.results);
figure(2)
plot_prediction(handles.in,handles.results);
figure(3)
plot_gaussian(handles.in,handles.results);
end

% --- Executes on button press in display_analysis.
function display_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to display_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

result_analysis(handles.in,handles.results);

if exist('log.txt','file')
    fid = fopen('log.txt');
    str = textscan(fid, '%s', 'Delimiter','\n'); str = str{1};
    fclose(fid);
    f=figure;
    hPan = uipanel(f,'Units','normalized');
    uicontrol(hPan, 'Style','listbox', ...
    'HorizontalAlignment','left', ...
    'Units','normalized', 'Position',[0 0 1 1], ...
    'String',str);
else
    set(handles.ErrorDisplay,'String','Could not file analysis log file');
    set(handles.ErrorDisplay,'ForegroundColor','r');
    set(hObject,'Enable','off')
end
    
end

% --- Executes on button press in showdata.
function showdata_Callback(hObject, eventdata, handles)
% hObject    handle to showdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure('position',[100 100 752 250]);
uitable(f,'unit','normalized','Position', [0 0 1 1],'Data',handles.in.data,'ColumnName',handles.in.headers);
end

% --- Executes on button press in feature_selection.
function feature_selection_Callback(hObject, eventdata, handles)
% hObject    handle to feature_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of feature_selection
end


% --- Executes on selection change in popup_xy.
function popup_xy_Callback(hObject, eventdata, handles)
% hObject    handle to popup_xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_xy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_xy

contents = cellstr(get(hObject,'String'));
index_sel = get(hObject,'Value');
value = contents(index_sel');
if strcmp(contents(index_sel'),'Most relevant')
    % do nothing
else
    % split string 
    strings = strsplit(value{1},' - ');
    handles.in.plot_xy=strings;
    handles.in.idx_xy = [];
    for i = 1:2
        idx = find(strcmp(handles.in.considered_inputs,handles.in.plot_xy{i}));
        if isempty(idx)
            error(['The specified axis input (' handles.in.plot_xy{i} ') could no be found in the list of selected input names']);
        end
        handles.in.idx_xy = [handles.in.idx_xy idx];
    end
    
end

guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function popup_xy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function ngrid_Callback(hObject, eventdata, handles)
% hObject    handle to ngrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.in.Ngrid = str2double(get(hObject,'String'));
set(handles.save_analysis,'Enable','off');

guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function ngrid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ngrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function description_Callback(hObject, eventdata, handles)
% hObject    handle to description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of description as text
%        str2double(get(hObject,'String')) returns contents of description as a double

handles.in.description = get(hObject,'String');

guidata(hObject,handles);

end

% --- Executes during object creation, after setting all properties.
function description_CreateFcn(hObject, eventdata, handles)
% hObject    handle to description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
