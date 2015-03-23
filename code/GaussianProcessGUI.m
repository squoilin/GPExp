function varargout = GaussianProcessGUI(varargin)
% GAUSSIANPROCESSGUI MATLAB code for GaussianProcessGUI.fig
%      GAUSSIANPROCESSGUI, by itself, creates a new GAUSSIANPROCESSGUI or raises the existing
%      singleton*.
%
%      H = GAUSSIANPROCESSGUI returns the handle to a new GAUSSIANPROCESSGUI or the handle to
%      the existing singleton*.
%
%      GAUSSIANPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAUSSIANPROCESSGUI.M with the given input arguments.
%
%      GAUSSIANPROCESSGUI('Property','Value',...) creates a new GAUSSIANPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GaussianProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GaussianProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GaussianProcessGUI

% Last Modified by GUIDE v2.5 03-Mar-2015 16:43:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GaussianProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GaussianProcessGUI_OutputFcn, ...
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


% --- Executes just before GaussianProcessGUI is made visible.
function GaussianProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GaussianProcessGUI (see VARARGIN)

% Choose default command line output for GaussianProcessGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GaussianProcessGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GaussianProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Select your file.
function SelectFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[str_filename,pathname,~] = uigetfile('*.*');

if ischar(str_filename) && ischar(pathname)
    
    error = 0;
    
    if strcmp(str_filename(end-2:end),'mat')
        load([pathname '/' str_filename]);
        vars = whos('-file',[pathname '/' str_filename]);
        
        if (sum(ismember({vars.name},'inputs'))~=1 || sum(ismember({vars.name},'outputs'))~=1 || sum(ismember({vars.name},'header_in'))~=1 || sum(ismember({vars.name},'header_out'))~=1)
            set(handles.ErrorDisplay,'String','The "in" structure should at least contain the "inputs" and "outputs" fields');

            ResetScript
            error = 1;
        else
            error = 0;
        end
        
    elseif strcmp(str_filename(end-2:end),'csv') || strcmp(str_filename(end-2:end),'xls') || strcmp(str_filename(end-3:end),'xlsx')
        %% XLS/CSV/XLSX import
        %The file must be organized into two sheets: one called 'inputs'
        %and the other called 'outputs'.
        %The header of each column will be considered as the variable name.
        %Each column represent then a variable and each line is one steady
        %state point.
        [~,sheets] = xlsfinfo([pathname '/' str_filename]);
        
        if sum(ismember(sheets,{'inputs' 'outputs'})) ~= 2
            set(handles.ErrorDisplay,'String','The file must contain the ''inputs'' and ''outputs'' sheets');
            ResetScript
            set(handles.ErrorDisplay,'ForegroundColor','r')
            error = 1;
        else
            set(handles.ErrorDisplay,'String','')
            
            [inputs,header_in,~] = xlsread([pathname '/' str_filename],'inputs');
            [outputs,header_out,~] = xlsread([pathname '/' str_filename],'outputs');
            if length(header_in)<1 || length(header_out)<1
                set(handles.ErrorDisplay,'String','Headers have not been provided.');
                ResetScript
                set(handles.ErrorDisplay,'ForegroundColor','r')
                error = 1;
            else
                error = 0;
            end
        end
    end
    
    if error == 0
        
        set(handles.ErrorDisplay,'String','')
        
        handles.in.inputs = inputs;
        handles.in.outputs = outputs;
        handles.n = size(handles.in.inputs,1);

        if size(handles.in.inputs,1) ~= size(handles.in.outputs,1)
            set(handles.ErrorDisplay,'String','The data matrix (inputs) and the targets (outputs) should have the same number of samples');
            ResetScript
            set(handles.ErrorDisplay,'ForegroundColor','r')
        else
            [m_in,n_in] = size(header_in);
            [m_out,n_out] = size(header_out);

            if m_in < n_in
                handles.in.inputnames = header_in';
            else
                handles.in.inputnames = header_in;
            end

            if m_out < n_out
                handles.in.outputnames = header_out';
            else
                handles.in.outputnames = header_out;
            end

            set(handles.ConsideredInputListbox,'Enable','on')
            set(handles.ConsideredInputListbox,'String',handles.in.inputnames)
            set(handles.ConsideredInputListbox,'Max',max(m_in,n_in))
            set(handles.ConsideredOutputListbox,'Enable','on')
            set(handles.ConsideredOutputListbox,'String',handles.in.outputnames)
            set(handles.ConsideredOutputListbox,'Max',1)%Only one single output can be selected at a time
            set(handles.TimeVariableCheckbox,'Enable','on')
            set(handles.kfoldTextbox,'Enable','on')
            set(handles.PermutationTextbox,'Enable','on')
            set(handles.RunButton,'Enable','on')
            set(handles.AnalysisStatusButton,'String','Run Analysis');
        end
    end
end

guidata(hObject,handles);

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

guidata(hObject,handles);

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

guidata(hObject,handles);

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


% --- Executes on button press in TimeVariableCheckbox.
function TimeVariableCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to TimeVariableCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TimeVariableCheckbox

handles.in.consider_temporality = get(hObject,'Value');

guidata(hObject,handles);



function kfoldTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to kfoldTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kfoldTextbox as text
%        str2double(get(hObject,'String')) returns contents of kfoldTextbox as a double

handles.in.kfolds = str2double(get(hObject,'String'));

% if ~isnan(str2double(get(handles.PermutationTextbox,'String')))
%     set(handles.RunButton,'Enable','on')
% end

guidata(hObject,handles);

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



function PermutationTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to PermutationTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PermutationTextbox as text
%        str2double(get(hObject,'String')) returns contents of PermutationTextbox as a double

handles.in.perm = str2double(get(hObject,'String'));

% if ~isnan(str2double(get(handles.kfoldTextbox,'String')))
%     set(handles.RunButton,'Enable','on')
% end

guidata(hObject,handles);

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


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.MainPlot)
cla(handles.SecondPlot)
cla(handles.ThirdPlot);
set( gcf, 'toolbar', 'figure' )

% handles.in = inputs_sanity_check(handles.in);

error = 0;

set(handles.ErrorDisplay,'String','');
set(handles.ErrorDisplay,'ForegroundColor','k');
set(handles.ResultTable,'Visible','off');
set(handles.Textbox1,'String','');
set(handles.Textbox2,'String','');

handles.in.kfolds = str2double(get(handles.kfoldTextbox,'String'));
handles.in.perm = str2double(get(handles.PermutationTextbox,'String'));

inputs_sanity_check_script;

if error == 0
    handles.results = main_model(handles.in);
    
    %% 3D plot of the results:
    x=handles.in.x;
    y=handles.in.y;
    z = handles.results.ndgrid.z;
    y_gp = handles.results.ndgrid.y_gp;
    s2 = handles.results.ndgrid.s2;

    nv = size(handles.in.considered_inputs,1);

    z_shaped = zeros(handles.in.Ngrid,nv);
    for i = 1:nv
        if nv > 1
            yp_shaped=reshape(y_gp,ones(1,nv)*handles.in.Ngrid);   % make a n-D matrix
            z_shaped(:,i)=linspace(min(z(:,i)),max(z(:,i)),handles.in.Ngrid);
        else
            z_shaped = z;
        end
    end

    axes(handles.MainPlot);

    if nv==1
        f = [y_gp+2*sqrt(s2); flipdim(y_gp-2*sqrt(s2),1)];
        fill([z; flipdim(z,1)], f, [7 7 7]/8)
        hold on; plot(z, y_gp); plot(x, y, '+');
        xlabel(handles.in.considered_inputs); ylabel(handles.in.considered_output);
        saveas(gcf,handles.in.name,'fig')
    elseif nv==2
        hold on;
        surf(z_shaped(:,1),z_shaped(:,2),yp_shaped);
        plot3(x(:,1),x(:,2), y, '+')
        xlabel(handles.in.considered_inputs(1)); ylabel(handles.in.considered_inputs(2)) ; zlabel(handles.in.considered_output);
        grid on
        view(45,25);
        saveas(gcf,handles.in.name,'fig')
    elseif nv > 2
        if isfield(handles.in,'idx_xy') %if the x and y axes are imposed
            idx_sort = [handles.in.idx_xy, 1:(min(handles.in.idx_xy)-1),(min(handles.in.idx_xy)+1):(max(handles.in.idx_xy)-1),(max(handles.in.idx_xy)+1:nv)] ;
        else % Sort the variables in terms of their weights or lengthscale in absolute value:
            if ~isempty(strfind(handles.in.covfunction{:},'ard'))
                [aa idx_sort] = sort(handles.results.hypcov(1:end-1));
            else
                [aa idx_sort] = sort(abs(handles.results.weights),'descend');
            end
        end
        % Set the nv-2 least relevant variables yp_shaped their median value:
        med = median(z);
        % Take a slice of the hypercube yp_shaped plot the two main relevant
        % variables:
        y_surf = permute(yp_shaped,idx_sort);
        med=med(idx_sort);
        itp=[];
        plottext = '';
        for i=1:nv
            if i==1 || i==2
                itp=[itp,',:'];
            else
                indi=idx_sort(i);
                vec = z_shaped(:,indi);
                [dd,indt] = min(abs(vec-med(i)));
                itp=[itp,[',',num2str(indt)]];
                plottext = strcat(plottext,{' '},handles.in.considered_inputs(indi), {' = '}, num2str(vec(indt)), {'; '});
            end
        end
        eval([ 'y_surf = y_surf(' itp(2:end) ');'])
        vecx = z_shaped(:,idx_sort(1));
        vecy = z_shaped(:,idx_sort(2));
        surf(vecy,vecx,y_surf);
        ylabel(handles.in.considered_inputs(idx_sort(1))); xlabel(handles.in.considered_inputs(idx_sort(2))) ; zlabel(handles.in.considered_output); %title(plottext);
        grid on
        view(45,25);
        saveas(gcf,handles.in.name,'fig')
    end

    %% plot of predicted vs. measured variables
    miny = min(y) - (max(y)-min(y))/15;
    maxy = max(y) + (max(y)-min(y))/15;

    axes(handles.SecondPlot);
    axis([miny,maxy,miny,maxy])
    hold on
    plot(y,handles.results.train.y_pred,'r*')
    plot([miny,maxy],[miny,maxy],'k')
    if isfield(handles.results,'CV')
        plot(y,handles.results.CV.y_pred,'b+')
        legend('Train (i.e. with all data samples)','45 deg','Cross-Validation','Location','best')
    else
        legend('Train (i.e. with all data samples)','45 deg','Location','best')
    end
    %title(['Predicted vs measured values of ' handles.in.considered_output])

    %% Plot all the points on the gaussian distribution + 5% confidence intervals

    axes(handles.ThirdPlot);
    alpha = 0.025;          % significance level
    mu = 0;               % mean
    sigma = 1;             % std
    range = max(3,max(abs(handles.results.outliers))) ;       % range (in the units of the std) for plotting
    cutoff1 = norminv(alpha, mu, sigma);
    cutoff2 = norminv(1-alpha, mu, sigma);
    x = [linspace(mu-range*sigma,cutoff1,20), ...
        linspace(cutoff1,cutoff2,50), ...
        linspace(cutoff2,mu+range*sigma,20)];
    y = normpdf(x, mu, sigma);
    plot(x,y)

    xlo = [-range x(x<=cutoff1) cutoff1];
    ylo = [0 y(x<=cutoff1) 0];
    patch(xlo, ylo, 'r','FaceAlpha',0.25)

    xhi = [cutoff2 x(x>=cutoff2) range];
    yhi = [0 y(x>=cutoff2) 0];
    patch(xhi, yhi, 'r','FaceAlpha',0.25)

    xlabel('Number of standard deviations')
    ylabel('Gaussian probability distribution')

    for i = 1:length(handles.results.outliers)
        string = num2str(i);
        xx = handles.results.outliers(i);
        yy=normpdf(xx,mu,sigma);
        text(xx,yy,string);
    end

    %title('Detection of outliers: Normal distribution of the error, with 5% significance intervals')

    
    set(handles.AnalysisStatusButton,'String','Analysis Done');
    
    %% Displaying the results
    set(handles.ResultTable,'Visible','on');
    if isfield(handles.results,'CV')
        set(handles.ResultTable,'ColumnName',{'Training Mode','Cross-Validation'});
        dataDisplayed = [handles.results.train.mae*100 handles.results.CV.mae*100;
                         handles.results.train.rsquare*100 handles.results.CV.rsquare*100;
                         handles.results.train.rmse handles.results.CV.rmse];
    else
        set(handles.ResultTable,'ColumnName',{'Training Mode'});
        dataDisplayed = [handles.results.train.mae*100;
                         handles.results.train.rsquare*100;
                         handles.results.train.rmse];
    end
    set(handles.ResultTable,'RowName',{'MARE [%]','R square [%]','RMSE'});
    set(handles.ResultTable,'Data',dataDisplayed);
    
    [val idx] = max(abs(handles.results.outliers)); 
    set(handles.Textbox1,'String',['Data point number ' num2str(idx) ', with an error ' num2str(val) ' times higher than the standard ' 'deviation of the gaussian process function at that particular point']);
    
    idx = find(abs(handles.results.outliers)>=1.96);
    text1=[];
    if isempty(idx)
        text1 = ['None'];
    else
        for i = 1:length(idx)
             text1=[text1 num2str(idx(i)) ', '];
        end
    end
    set(handles.Textbox2,'String',{['The following data points present a significance level lower than 5 percent.' 'They are therefore likely to be outliers:'], ['Data point number ' text1(1:end-1)]});
    
end

% --- Executes during object creation, after setting all properties.
function MainPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MainPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate MainPlot


% --- Executes when entered data in editable cell(s) in ResultTable.
function ResultTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ResultTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in ResultTable.
function ResultTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to ResultTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AutoRunButton.
function AutoRunButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutoRunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N_var = length(handles.in.inputnames);
N_comb = factorial(N_var)./(factorial([N_var-1:-1:0]).*factorial([1:N_var])); %Number of possible cases per level (ie for 1 input, 2 inputs, 3 inputs ...)
N_comb_tot = sum(N_comb);

for i = 1%:N_var %Number of variable "trees"
    for j = 1:N_var-i
        
    end
end

