error_check = 0;

while error_check~=1
    if ~isfield(handles.in,'inputnames') %name of the input variables
        handles.in.inputnames=cell(size(in.inputs,2),1);
        for i=1:size(handles.in.inputs,2)
            handles.in.inputnames{i} = ['Input ',num2str(i)];
        end
    else
        if ~isvector(handles.in.inputnames)
            set(handles.ErrorDisplay,'String','The inputnames variable must be a vector');
            set(handles.ErrorDisplay,'ForegroundColor','r');
            error = 1;
            break
        end
        if length(handles.in.inputnames) ~= size(handles.in.inputs,2)
            set(handles.ErrorDisplay,'String',['The number of specified input names (' num2str(size(handles.in.inputnames,1)) ') does not match with the number of columns of the input array (' num2str(size(handles.in.inputs,2)) ')']);
            set(handles.ErrorDisplay,'ForegroundColor','r');
            error = 1;
            break
        end
    end

    if ~isfield(handles.in,'outputnames') %name of the output variables
        handles.in.outputnames=cell(size(handles.in.outputs,2),1);
        for i=1:size(handles.in.outputs,2)
            handles.in.outputnames{i} = ['output ',num2str(i)];
        end
    else
        if size(handles.in.outputnames,1) ~= size(handles.in.outputs,2)
            set(handles.ErrorDisplay,'String',['The number of specified output names (' num2str(size(handles.in.outputnames,1)) ') does not match with the number of columns of the output array (' num2str(size(handles.in.outputs,2)) ')']);
            set(handles.ErrorDisplay,'ForegroundColor','r');
            error = 1;
            break
        end
    end

    list_inputs = [];
    for i=1:size(handles.in.inputs,2)
        list_inputs = [list_inputs char(10) handles.in.inputnames{i}];
    end

    list_outputs = [];
    for i=1:size(handles.in.outputs,2)
        list_outputs = [list_outputs char(10) handles.in.outputnames{i}];
    end

    % Checking output and output name:
    if ~isfield(handles.in,'considered_output') && size(handles.in.outputs,2) ~= 1
        set(handles.ErrorDisplay,'String',['Please specify the output to be considered in the "considered_output" variable. Possible values are: ' list_outputs]);
        set(handles.ErrorDisplay,'ForegroundColor','r');
        error = 1;
        break
    elseif ~isfield(handles.in,'considered_output') && size(handles.in.output,2) == 1
        handles.in.considered_output = handles.in.outputnames{1};
    else
        if iscell(handles.in.considered_output) && size(handles.in.considered_output,1)>1
            set(handles.ErrorDisplay,'String','Predicting one target (output) at a time, please execute multiple times for multiple targets');
            set(handles.ErrorDisplay,'ForegroundColor','r');
            error = 1;
            break
        end
        idx_output = find(strcmp(handles.in.outputnames,handles.in.considered_output));
        if isempty(idx_output)
            set(handles.ErrorDisplay,'String',['The specified output (' handles.in.considered_output ') could no be found in the list of output names:' list_outputs]);
            set(handles.ErrorDisplay,'ForegroundColor','r');
            error = 1;
            break
        end
        handles.in.y = handles.in.outputs(:,idx_output);
    end
    
    % Checking inputs and inputs names:
    if ~isfield(handles.in,'considered_inputs') && size(handles.in.inputs,2) ~= 1
        set(handles.ErrorDisplay,'String',['Please specify the inputs to be considered in the "considered_inputs" variable. Possible values are: ' list_inputs]);
        set(handles.ErrorDisplay,'ForegroundColor','r');
        error = 1;
        break
    elseif ~isfield(handles.in,'considered_inputs') && size(handles.in.input,2) == 1
        handles.in.considered_inputs = handles.in.inputnames{1};
    else
        idx_inputs = [];
        for i = 1:size(handles.in.considered_inputs,1)
            idx_input = find(strcmp(handles.in.inputnames,handles.in.considered_inputs{i}));
            if isempty(idx_input)
                set(handles.ErrorDisplay,'String',['The specified input (' handles.in.considered_inputs{i} ') could no be found in the list of input names:' list_inputs]);
                set(handles.ErrorDisplay,'ForegroundColor','r');
                error = 1;
                break
            end
            idx_inputs = [idx_inputs idx_input];
        end
        handles.in.x = handles.in.inputs(:,idx_inputs);
    end

    if ~isfield(handles.in,'consider_temporality')
        handles.in.consider_temporality = false;
    end

    if handles.in.consider_temporality && sum(strcmp('idx_{sample}',handles.in.considered_inputs))<1
        handles.in.x = horzcat(handles.in.x,(1:size(handles.in.x,1))');
        handles.in.considered_inputs{end+1,1}='idx_{sample}';
    end
    
    Nrows = size(handles.in.x,1);
    Nvars = size(handles.in.x,2);

    list_considered_inputs = [];
    for i=1:size(handles.in.considered_inputs,1)
        list_considered_inputs = [list_considered_inputs char(10) handles.in.considered_inputs{i}];
    end

    if isfield(handles.in,'plot_xy')
        if size(handles.in.plot_xy,1) ~= 2
            set(handles.ErrorDisplay,'String','plot_xy should contain exactly two cells (one for the x axis and one for the y axis)');
            set(handles.ErrorDisplay,'ForegroundColor','r');
            error = 1;
            break
        else
            handles.in.idx_xy = [];
            for i = 1:2
                idx = find(strcmp(handles.in.considered_inputs,handles.in.plot_xy{i}));
                if isempty(idx)
                    set(handles.ErrorDisplay,'String',['The specified axis input (' handles.in.plot_xy{i} ') could no be found in the list of selected input names:' list_considered_inputs]);
                    set(handles.ErrorDisplay,'ForegroundColor','r');
                    error = 1;
                    break
                end
                handles.in.idx_xy = [handles.in.idx_xy idx];
            end
        end
    end

    if ~isfield(handles.in,'covfunction') %covariance function
        handles.in.covfunction={'covSEard'};
    end

    if ~isfield(handles.in,'kfolds') %number of folds for cross-validation
        in.kfolds=0;
    end

    if ~isfield(handles.in,'perm') %number of permutations for significance
        handles.in.perm=0;
    end
    disp(handles.in.kfolds)
    if handles.in.kfolds < 0 && handles.in.perm > 0
        set(handles.ErrorDisplay,'String','Permutations only work with cross validation. Please set the number of folds to a non-negative value or set the number of permutations to zero');
        set(handles.ErrorDisplay,'ForegroundColor','r');
        error = 1;
        break
    end

    if ~isfield(handles.in,'name') %name of the simulation for plotting and for saving results
        handles.in.name = 'default';
    end
    
    if isfield(handles.in,'description')
        if ~iscellstr(handles.in.description)
            set(handles.ErrorDisplay,'String','The "description" field should be a cell array of strings');
            set(handles.ErrorDisplay,'ForegroundColor','r');
            error = 1;
            break
        end
    else
        handles.in.description = {' '};
    end

    if ~isfield(handles.in,'Ngrid') %number of grid points per dimension
        Ngrid_max = floor(1E5^(1/Nvars));
        handles.in.Ngrid = min(20,Ngrid_max);
    end

    Ngridpoints = handles.in.Ngrid^Nvars;
    if Ngridpoints > 1E6
        set(handles.ErrorDisplay,'String','You have many grid points (probably due to a high number of considered input variables). This might cause memory issues and can be solved by reducing the Ngrid parameter')
    end
    
    if ~isfield(handles.in,'hyp')
        handles.in.hyp.cov = zeros(Nvars+1,1);
        handles.in.hyp.lik = 0;
    end

    handles.in.checked = true;
    
    error = 0;
    error_check = 1;
end